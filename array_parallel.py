import os, fcntl, sys
import numpy as np
from scipy import linalg as sla

# --- CLI / SLURM ---
task_id = int(sys.argv[1])                       # array index
total_tasks = int(os.environ.get("SLURM_ARRAY_TASK_COUNT", "1"))
job_id = os.environ.get("SLURM_JOB_ID", "local")

# --- model params ---
cell_size = 3           # 3 or 4
extra_sites = True
N_cells = 50
Nt = N_cells * (cell_size**2)
v = 0.2
delta = 0.4 * np.pi
alpha = (np.sqrt(5.0) - 1.0) / 2.0
t = 1.0
lm = 1                  # nn hop

# --- helper: choose energy window on Re(E) ---
def energy_window(sel):
    if isinstance(sel, (tuple, list)) and len(sel) == 2:
        lo, hi = float(sel[0]), float(sel[1]); return lo, hi
    if sel in ("3x3", 3): return -2.0, -1.0
    if sel in ("4x4", 4): return -3.0, -2.0
    return -3.0, -2.0

# --- coordinates ---
def make_coordinates(Nt, N_):
    custom_rangey = np.arange(start=1 - (N_ + 1) / 2, stop=(N_ - 1) / 2 + 0.01)
    if extra_sites:
        if N_ % 2 != 0:
            Y_base = np.array([0.0] + [y for _ in range(N_) for y in custom_rangey])
            Y = [y for _ in range(N_cells) for y in Y_base] + [0.0]
        else:
            Y_base = np.array([-0.5] + [0.5] + [y for _ in range(N_) for y in custom_rangey])
            Y = [y for _ in range(N_cells) for y in Y_base] + [-0.5, 0.5]
    else:
        Y_base = np.array([y for _ in range(N_) for y in custom_rangey])
        Y = [y for _ in range(N_cells) for y in Y_base]

    if extra_sites:
        rng = np.linspace(-N_cells * (cell_size + 1) / 2.0,
                          N_cells * (cell_size + 1) / 2.0,
                          N_cells * cell_size + N_cells + 1)
        X = []
        lk = 1 if (N_ % 2 != 0) else 2
        for i in range(len(rng)):
            if i % (cell_size + 1) == 0:
                X.extend([rng[i]] * lk)
            else:
                X.extend([rng[i]] * cell_size)
    else:
        rng = np.linspace(-N_cells * (cell_size - 1) / 2.0,
                          N_cells * (cell_size - 1) / 2.0,
                          N_cells * cell_size + N_cells - 1)
        X = []
        lk = 1 if (N_ % 2 != 0) else 2
        for i in range(len(rng)):
            X.extend([rng[i]] * cell_size)

    X = np.array(X, dtype=float)
    Y = np.array(Y, dtype=float)
    Nl = len(Y)
    return X, Y, custom_rangey, Nl

X, Y, custom_range, Nl = make_coordinates(Nt, cell_size)

# --- phases / algebra helpers ---
def phase_factor(x1, y1, x2, y2, B):
    return np.exp(-1j * 2.0 * np.pi * B * (x2 - x1) * ((y1 + y2) / 2.0))

def orthomatrix(vl1, vl2, Nl):
    M = np.zeros((Nl, Nl), dtype=np.complex128)
    for i in range(Nl):
        for j in range(Nl):
            M[i, j] = np.dot(np.conj(vl1[:, i]), vl2[:, j])
    return M

def newbiorth1(vl1, vl2):
    M = orthomatrix(vl1, vl2, Nl)
    P, L, U = sla.lu(M)
    L1 = P @ L
    vlp = sla.inv(L1) @ np.conj(vl1.T)
    vrp = vl2 @ sla.inv(U)
    return np.conj(vlp.T), vrp

def ipr(vl1, vl2):
    overlaps = np.conj(vl1) * vl2
    return np.sum(np.abs(overlaps) ** 2, axis=0)

# --- build H (non-Hermitian via imaginary onsite) ---
def eigsys(X,Y,v, N_, t, B, Nl):
    H = np.zeros((Nl, Nl), dtype=np.complex128)

    for i, (x, y) in enumerate(zip(X, Y)):
        for j, (x1, y1) in enumerate(zip(X, Y)):
            # vertical
            if (x1 == x) and (y1 == y + lm):
                val = t * phase_factor(x, y, x1, y1, B)
                H[i, j] = val
                H[j, i] = np.conj(val)
            # horizontal
            if (y1 == y) and (x1 == x + lm):
                val = t * phase_factor(x, y, x1, y1, B)
                H[i, j] = val
                H[j, i] = np.conj(val)

    if extra_sites and (N_ % 2 == 0):
        for k0 in range(0, Nl, cell_size**2 + 2):
            H[k0, k0 + 1] = 0.0
            H[k0 + 1, k0] = 0.0

    # imaginary onsite → non-Hermitian
    k = 0
    if extra_sites:
        if N_ % 2 == 0:
            for i0 in range(2, len(X) - 1, cell_size**2 + 2):
                for j in range(i0, i0 + cell_size**2):
                    H[j, j] = 1j * v * np.sin(2 * np.pi * alpha * k + delta)
                k += 1
        else:
            for i0 in range(1, len(X) - 1, cell_size**2 + 1):
                for j in range(i0, i0 + cell_size**2):
                    H[j, j] = 1j * v * np.sin(2 * np.pi * alpha * k + delta)
                k += 1
    return H

# --- metric ---
def compute_percentage(B, v, cell_size, t, Nl, window=cell_size):
    H = eigsys(X,Y,v, cell_size, t, B, Nl)
    # non-Hermitian eigensystem
    w, VL, VR = sla.eig(H, left=True, right=True)
    idx = np.argsort(w.real)
    w, VL, VR = w[idx], VL[:, idx], VR[:, idx]

    VLp, VRp = newbiorth1(VL, VR)

    e_lo, e_hi = energy_window(window)
    mask = (w.real >= e_lo) & (w.real <= e_hi)
    total = int(np.sum(mask))
    if total == 0:
        return (B, v, 0.0)

    f = ipr(VLp[:, mask], VRp[:, mask])
    thr = 12.0 / Nl
    percentage = 100.0 * np.sum(np.abs(f) < thr) / total
    return (B, v, float(percentage))

# --- parameter grids ---
B_vals = np.linspace(0.0025, 0.025, 100)
v_vals = np.linspace(0.0, 0.2,   100)

def run_block_for_task(task_id, total_tasks, B_vals, v_vals, cell_size, t, Nl):
    rows = []
    if total_tasks == 1:
        iterable = ((B, v) for B in B_vals for v in v_vals)
    else:
        iterable = []
        if 0 <= task_id < len(B_vals):
            B = B_vals[task_id]
            iterable = ((B, v) for v in v_vals)
    for B, v in iterable:
        rows.append(compute_percentage(B, v, cell_size, t, Nl, window=cell_size))
    return rows

def append_rows_atomic(path, rows, header="B v percentage"):
    if not rows:
        return
    A = np.asarray(rows, dtype=float)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "a+") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            f.seek(0, os.SEEK_END)
            if f.tell() == 0:
                f.write(header + "\n")
            np.savetxt(f, A, fmt="%.8g")
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)

# --- run & append ---
rows = run_block_for_task(task_id, total_tasks, B_vals, v_vals, cell_size, t, Nl)
out_file = f"results/imag_even_parallel_array_ALL_{job_id}.dat"
append_rows_atomic(out_file, rows)
print(f"Task {task_id}: appended {len(rows)} rows to {out_file}")

