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

