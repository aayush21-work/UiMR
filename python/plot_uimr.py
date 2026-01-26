import numpy as np
import matplotlib.pyplot as plt

def visualize_macro_blocks(filename, n_blocks_row, n_blocks_col):
    """
    Visualizes the bounds.dat file where each line contains: 
    le ri to bo
    """
    try:
        # Load the data: Each row is [le, ri, to, bo]
        data = np.loadtxt(filename)
    except OSError:
        print(f"Error: Could not find {filename}.")
        return

    # Check data integrity
    expected_lines = n_blocks_row * n_blocks_col
    if data.shape[0] != expected_lines:
        print(f"Data Mismatch! Expected {expected_lines} blocks, found {data.shape[0]}.")
        print("Check your 'rows' and 'cols' variables in the Python script.")
        return

    # 1. Reshape data to represent the 2D grid of blocks
    # Shape becomes: (Rows, Cols, 4_values)
    grid_data = data.reshape((n_blocks_row, n_blocks_col, 4))

    # 2. Calculate the average potential of each block
    # We simply average (Left + Right + Top + Bottom) / 4 to get a 'center' value estimate
    block_averages = np.mean(grid_data, axis=2)

    # 3. Setup the Plot
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot the heatmap
    im = ax.imshow(block_averages, cmap='inferno', origin='upper')
    
    # Add colorbar
    plt.colorbar(im, label='Average Potential per Block')

    # 4. Annotate the grid with the values (Optional, good for debugging)
    # We will print the value inside each square
    for i in range(n_blocks_row):
        for j in range(n_blocks_col):
            val = block_averages[i, j]
            # Choose text color based on brightness for readability
            text_color = "black" if val > np.max(block_averages)/2 else "white"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center", color=text_color, fontsize=8)

    # Labels and Titles
    ax.set_title(f"Macro-Grid Potentials\n(Grid Size: {n_blocks_row}x{n_blocks_col})")
    ax.set_xlabel("Block Column Index")
    ax.set_ylabel("Block Row Index")
    
    # Ticks setup to look like a grid
    ax.set_xticks(np.arange(n_blocks_col))
    ax.set_yticks(np.arange(n_blocks_row))
    ax.set_xticklabels(np.arange(1, n_blocks_col+1))
    ax.set_yticklabels(np.arange(1, n_blocks_row+1))

    plt.tight_layout()
    plt.show()

# ==========================================
# CONFIGURATION
# ==========================================
# Based on your uimr(2, ...) call where rows=10, cols=10:
# rows becomes 10/2 = 5
# cols becomes 10/2 = 5
# So you have a 5x5 grid of blocks.

BLOCK_ROWS = 5 
BLOCK_COLS = 5

visualize_macro_blocks("bounds.dat", BLOCK_ROWS, BLOCK_COLS)
