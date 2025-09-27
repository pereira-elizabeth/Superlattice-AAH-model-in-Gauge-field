def compute_percentage(B, v, cell_size, t, Nl, window):
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
