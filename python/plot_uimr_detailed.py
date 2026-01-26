import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

def visualize_detailed_boundaries(filename, n_blocks_row, n_blocks_col):
    """
    Visualizes bounds.dat by splitting each block into 4 triangles:
    Left, Right, Top, Bottom.
    """
    try:
        # Load data. Format: le, ri, to, bo
        raw_data = np.loadtxt(filename)
    except OSError:
        print(f"Error: Could not find {filename}.")
        return

    # reshape to (Rows, Cols, 4)
    if raw_data.shape[0] != n_blocks_row * n_blocks_col:
        print(f"Mismatch: Expected {n_blocks_row*n_blocks_col} rows, got {raw_data.shape[0]}")
        return
        
    grid_data = raw_data.reshape((n_blocks_row, n_blocks_col, 4))

    # Setup Plot
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # We will collect all triangles here to plot them efficiently
    triangle_patches = []
    triangle_values = []

    # Loop through every block in the grid
    for r in range(n_blocks_row):
        for c in range(n_blocks_col):
            # Extract values for this block
            # Order in C++ file: le (0), ri (1), to (2), bo (3)
            val_l = grid_data[r, c, 0]
            val_r = grid_data[r, c, 1]
            val_t = grid_data[r, c, 2]
            val_b = grid_data[r, c, 3]

            # Coordinates for the corners of the block
            # (Using Matrix coordinates: Y grows downwards)
            x_left = c
            x_right = c + 1
            y_top = r
            y_bottom = r + 1
            
            x_center = c + 0.5
            y_center = r + 0.5

            # Define the 4 Triangles (Polygon vertices)
            
            # 1. Left Triangle (Vertices: Top-Left, Bottom-Left, Center)
            tri_l = [[x_left, y_top], [x_left, y_bottom], [x_center, y_center]]
            
            # 2. Right Triangle (Vertices: Top-Right, Bottom-Right, Center)
            tri_r = [[x_right, y_top], [x_right, y_bottom], [x_center, y_center]]
            
            # 3. Top Triangle (Vertices: Top-Left, Top-Right, Center)
            tri_t = [[x_left, y_top], [x_right, y_top], [x_center, y_center]]
            
            # 4. Bottom Triangle (Vertices: Bottom-Left, Bottom-Right, Center)
            tri_b = [[x_left, y_bottom], [x_right, y_bottom], [x_center, y_center]]

            # Add to collection
            triangle_patches.extend([
                patches.Polygon(tri_l), 
                patches.Polygon(tri_r), 
                patches.Polygon(tri_t), 
                patches.Polygon(tri_b)
            ])
            
            # Store values in same order
            triangle_values.extend([val_l, val_r, val_t, val_b])

    # Create the collection
    p = PatchCollection(triangle_patches, cmap='inferno', alpha=1.0)
    p.set_array(np.array(triangle_values))
    p.set_edgecolor('face') # Removes lines between triangles for smoother look
    
    ax.add_collection(p)
    
    # Add colorbar
    plt.colorbar(p, label='Boundary Potential')

    # Formatting the axes to look like a matrix (0,0 at top-left)
    ax.set_xlim(0, n_blocks_col)
    ax.set_ylim(n_blocks_row, 0) # Invert Y to match matrix index
    ax.set_aspect('equal')
    
    # Add grid lines to separate blocks clearly
    ax.set_xticks(np.arange(n_blocks_col + 1))
    ax.set_yticks(np.arange(n_blocks_row + 1))
    ax.grid(color='white', linestyle='-', linewidth=0.5, alpha=0.5)

    ax.set_title(f"Detailed Boundary Visualization\n(Left/Right/Top/Bottom per Block)")
    plt.show()

# ============================
# CONFIGURATION
# ============================
BLOCK_ROWS = 5
BLOCK_COLS = 5

visualize_detailed_boundaries("bounds.dat", BLOCK_ROWS, BLOCK_COLS)
