import numpy as np
from scipy import linalg as sla

def energy_window(sel):
    if isinstance(sel, (tuple, list)) and len(sel) == 2:
        lo, hi = float(sel[0]), float(sel[1]); return lo, hi
    if sel in ("3x3", 3): return -2.0, -1.0
    if sel in ("4x4", 4): return -3.0, -2.0
    return -3.0, -2.0

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
def eigsys(X,Y,v, N_, t, B, Nl, lm):
    extra_sites = True
    cell_size = N_
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
