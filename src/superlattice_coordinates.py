import numpy as np

def make_coordinates(Nt, N_, N_cells, extra_sites = True):
    cell_size = N_
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

