import numpy as np
import matplotlib.pyplot as plt
import glob

# ================= USER SETTINGS =================
Nb_x = 10      # blocks in x
Nb_y = 10      # blocks in y

rows = 10       # rows per block (same as C++)
cols = 10     # cols per block (same as C++)

pattern = "sol_*.dat"
# =================================================


# collect and sort solution files
files = sorted(
    glob.glob(pattern),
    key=lambda s: int(s.split("_")[1].split(".")[0])
)

assert len(files) == Nb_x * Nb_y, \
    f"Expected {Nb_x*Nb_y} blocks, found {len(files)}"

# allocate global field
G = np.zeros((Nb_y * rows, Nb_x * cols))

k = 0
for by in range(Nb_y):
    for bx in range(Nb_x):

        fname = files[k]
        B = np.loadtxt(fname)

        assert B.shape == (rows, cols), \
            f"Shape mismatch in {fname}"

        i0 = by * rows
        j0 = bx * cols

        G[i0:i0+rows, j0:j0+cols] = B
        k += 1


# ===================== PLOT ======================
plt.figure(figsize=(7, 6))

plt.imshow(G, origin="lower", cmap="inferno")
plt.colorbar(label="Potential")

plt.title("Stitched Laplace Solution")
plt.xlabel("x")
plt.ylabel("y")

# draw block grid (optional but useful)
for y in range(0, G.shape[0], rows):
    plt.axhline(y - 0.5, color="white", linewidth=0.4)

for x in range(0, G.shape[1], cols):
    plt.axvline(x - 0.5, color="white", linewidth=0.4)

plt.tight_layout()
plt.show()
