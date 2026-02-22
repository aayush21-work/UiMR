import numpy as np
import matplotlib.pyplot as plt
import glob, os

TOTAL_ROWS = 100
TOTAL_COLS = 100
N          = 10         # change as per the cpp or the cuda file
sub_rows   = TOTAL_ROWS // N
sub_cols   = TOTAL_COLS // N

full = np.zeros((TOTAL_ROWS, TOTAL_COLS))

n_files = N * N
for sid in range(n_files):
    fname = f"sol_{sid}.dat"
    if not os.path.exists(fname):
        raise FileNotFoundError(f"Missing {fname} — make sure you run the C++ first.")
    
    sub = np.loadtxt(fname)        
    # subdomain id maps to tile position
    # id increments j (col) fastest, then i (row) — matching uimr_init loop
    tile_i = sid // N               # row-tile index
    tile_j = sid %  N               # col-tile index

    r0 = tile_i * sub_rows
    c0 = tile_j * sub_cols
    full[r0:r0+sub_rows, c0:c0+sub_cols] = sub


fig, ax = plt.subplots(figsize=(7, 6))

im = ax.imshow(full, origin="upper", cmap="hot",
               extent=[0, TOTAL_COLS, TOTAL_ROWS, 0])

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Tempreture", fontsize=12)

ax.set_title(f"Laplace Solution — {N}×{N} domain decomposition", fontsize=13)
ax.set_xlabel("x")
ax.set_ylabel("y")

for k in range(1, N):
    ax.axhline(k * sub_rows, color="cyan", lw=0.4, alpha=0.5)
    ax.axvline(k * sub_cols, color="cyan", lw=0.4, alpha=0.5)

plt.tight_layout()
plt.show()
