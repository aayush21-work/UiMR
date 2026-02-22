import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# ── parameters — must match your C++ solver ──────────────────────────────────
ROWS = 500
COLS = 500
N    = 25          # grid of N×N tiles
SR   = ROWS // N   # tile rows
SC   = COLS // N   # tile cols
NSUB = N * N
# ─────────────────────────────────────────────────────────────────────────────

def load_uimr_omp(filename="../output/sol_uimr_omp.dat"):
    """Load the single merged OMP output file."""
    data = np.loadtxt(filename)   # shape: (ROWS, COLS)
    assert data.shape == (ROWS, COLS), \
        f"Expected ({ROWS},{COLS}), got {data.shape}"
    return data

def stitch_tiles(prefix="sol_", nsub=NSUB, n=N, sr=SR, sc=SC):
    """Stitch N×N individual sol_<id>.dat tile files into full domain."""
    full = np.zeros((ROWS, COLS))
    for tid in range(nsub):
        fname = f"{prefix}{tid}.dat"
        if not os.path.exists(fname):
            raise FileNotFoundError(f"Missing tile file: {fname}")
        tile = np.loadtxt(fname)   # shape: (sr, sc)
        ti = tid // n              # row index of tile in grid
        tj = tid  % n              # col index of tile in grid
        full[ti*sr:(ti+1)*sr, tj*sc:(tj+1)*sc] = tile
    return full

def load_traditional(filename="sol_traditional.dat"):
    """Load the traditional single-file solution."""
    data = np.loadtxt(filename)
    assert data.shape == (ROWS, COLS), \
        f"Expected ({ROWS},{COLS}), got {data.shape}"
    return data

def plot_solution(data, title, ax, vmin=None, vmax=None):
    im = ax.imshow(data, origin="upper", cmap="hot",
                   vmin=vmin, vmax=vmax, aspect="auto",
                   extent=[0, COLS, ROWS, 0])
    ax.set_title(title, fontsize=13)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return im

# ── load OMP merged file only ─────────────────────────────────────────────────
if not os.path.exists("sol_uimr_omp.dat"):
    print("sol_uimr_omp.dat not found. Run the OMP solver first.")
    sys.exit(1)

data = load_uimr_omp("sol_uimr_omp.dat")
print("Loaded sol_uimr_omp.dat")

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(data, origin="upper", cmap="hot", aspect="auto",
               extent=[0, COLS, ROWS, 0])
ax.set_title("Laplace Solution — UIMR OMP", fontsize=14)
ax.set_xlabel("x")
ax.set_ylabel("y")
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Temperature", fontsize=11)

fig.tight_layout()
fig.savefig("solution_plot.png", dpi=150)
print("Saved solution_plot.png")
plt.show()
