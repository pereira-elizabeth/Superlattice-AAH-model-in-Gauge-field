# Make repo root importable (works in CI and locally)
import sys, pathlib, os, tempfile
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np

from src.superlattice_coordinates import make_coordinates
from src.parallelsource import run_block_for_task, append_rows_atomic

# If your module name is different, fix these imports accordingly:
from src.hamiltonian_gauge import eigsys, ipr, newbiorth1, energy_window, phase_factor, orthomatrix
from src.percentage_metric import compute_percentage


def test_make_coordinates_tiny_grid():
    """Coordinates function returns expected shapes for a tiny system."""
    cell_size = 3
    N_cells = 2                      # tiny, fast
    Nt = N_cells * (cell_size**2)    # must be consistent with your code

    X, Y, custom_range, Nl = make_coordinates(Nt, cell_size, N_cells)

    # Basic sanity checks
    assert isinstance(Nl, int) and Nl > 0
    assert np.asarray(X).size == Nt
    assert np.asarray(Y).size == Nt
    assert hasattr(custom_range, "__iter__")  # any iterable is fine


def test_run_block_for_task_returns_rows():
    """
    run_block_for_task should accept small grids and return rows like [B, v, percentage].
    We keep the grid very small so CI stays fast.
    """
    cell_size = 3
    N_cells = 2
    Nt = N_cells * (cell_size**2)

    # We only need Nl for run_block_for_task; get it via the same helper the script uses.
    _, _, _, Nl = make_coordinates(Nt, cell_size, N_cells)

    # Tiny parameter sets
    B_vals = np.array([0.005, 0.010])
    v_vals = np.array([0.0, 0.1])

    # Single task partition (no SLURM)
    task_id = 0
    total_tasks = 1

    rows = run_block_for_task(task_id, total_tasks, B_vals, v_vals, cell_size, 1.0, Nl)

    # rows should be an iterable of tuples/lists; each row should be length >= 3
    rows = list(rows)
    assert len(rows) >= 1
    for r in rows:
        assert len(r) >= 3
        B, v, pct = float(r[0]), float(r[1]), float(r[2])
        assert np.isfinite(B) and np.isfinite(v) and np.isfinite(pct)
        # percentage should be in a sane range
        assert -1e-6 <= pct <= 100 + 1e-6


def test_append_rows_atomic_writes_all_rows(tmp_path=None):
    """
    append_rows_atomic should create/append lines atomically.
    We write a few rows to a temp file and confirm line count matches.
    """
    cell_size = 3
    N_cells = 2
    Nt = N_cells * (cell_size**2)
    _, _, _, Nl = make_coordinates(Nt, cell_size, N_cells)

    B_vals = np.array([0.005, 0.010])
    v_vals = np.array([0.0, 0.1])
    rows = list(run_block_for_task(0, 1, B_vals, v_vals, cell_size, 1.0, Nl))

    # Temp output file under tests/temp
    if tmp_path is None:
        tmpdir = tempfile.TemporaryDirectory()
        out_file = os.path.join(tmpdir.name, "test_output.dat")
    else:
        out_file = os.path.join(tmp_path, "test_output.dat")

    # First append
    append_rows_atomic(out_file, rows)
    # Append again to ensure it appends (not overwrite)
    append_rows_atomic(out_file, rows)

    with open(out_file, "r") as f:
        lines = [ln for ln in f.read().splitlines() if ln.strip()]

    assert len(lines) == 2 * len(rows)
