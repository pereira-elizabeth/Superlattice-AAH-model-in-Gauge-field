import os, fcntl, sys
import numpy as np
from scipy import linalg as sla
from src.hamiltonian_gauge import energy_window, phase_factor, orthomatrix, newbiorth1, ipr, eigsys
from src.superlattice_coordinates import make_coordinates
from src.percentage_matric import compute_percentage
from src.parallelsource import append_rows_atomic, run_block_for_task

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

X, Y, custom_range, Nl = make_coordinates(Nt, cell_size)

# --- parameter grids ---
B_vals = np.linspace(0.0025, 0.025, 100)
v_vals = np.linspace(0.0, 0.2,   100)

# --- run & append ---
rows = run_block_for_task(task_id, total_tasks, B_vals, v_vals, cell_size, t, Nl)
out_file = f"results/imag_even_parallel_array_ALL_{job_id}.dat"
append_rows_atomic(out_file, rows)
print(f"Task {task_id}: appended {len(rows)} rows to {out_file}")

