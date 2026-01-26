import numpy as np
import matplotlib.pyplot as plt

# ================= USER SETTINGS =================
bounds_file = "bounds.dat"

Nb_x = 5     # number of blocks in x
Nb_y = 5     # number of blocks in y

rows = 40    # visual resolution per block (height)
cols = 40    # visual resolution per block (width)
# =================================================


# load bounds: each row = [left, right, top, bottom]
bounds = np.loadtxt(bounds_file)

assert bounds.shape[0] == Nb_x * Nb_y, \
    "Mismatch between bounds.dat and block layout"

# global canvas (NaN interior)
G = np.full((Nb_y * rows, Nb_x * cols), np.nan)

k = 0
for by in range(Nb_y):
    for bx in range(Nb_x):

        le, ri, to, bo = bounds[k]
        k += 1

        i0 = by * rows
        j0 = bx * cols

        # paint boundaries
        G[i0:i0+rows, j0]           = le    # left
        G[i0:i0+rows, j0+cols-1]    = ri    # right
        G[i0, j0:j0+cols]           = bo    # bottom
        G[i0+rows-1, j0:j0+cols]    = to    # top


# ===================== PLOT ======================
plt.figure(figsize=(8, 6))

cmap = plt.cm.inferno.copy()
cmap.set_bad(color="white")  # NaN interior

plt.imshow(G, origin="lower", cmap=cmap)
plt.colorbar(label="Boundary value")

plt.title("Stitched Block Boundaries")
plt.xlabel("x")
plt.ylabel("y")

# block grid
for y in range(0, G.shape[0], rows):
    plt.axhline(y - 0.5, color="black", linewidth=0.4)

for x in range(0, G.shape[1], cols):
    plt.axvline(x - 0.5, color="black", linewidth=0.4)

plt.tight_layout()
plt.show()
