#include <iostream>
#include <fstream>
#include <vector>

constexpr int rows = 5000;
constexpr int cols = 5000;
constexpr int nt   = 5000;

/* ---------- index helper ---------- */
inline int idx(int k, int i, int j) {
    return k * rows * cols + i * cols + j;
}

/* ---------- write history ---------- */
void write_history(const std::vector<double>& history) {
    std::ofstream file("solution.dat");

    for (int k = 0; k < nt; k++) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                file << history[idx(k, i, j)] << " ";
            }
            file << "\n";
        }
        file << "\n"; // separate time slices
    }
}

int main() {
    /* ---------- main field ---------- */
    double A[rows][cols] = {};

    /* ---------- history on heap ---------- */
    std::vector<double> history(nt * rows * cols, 0.0);

    /* ---------- boundary conditions ---------- */
    for (int i = 0; i < rows; i++) {
        A[i][0]        =100.0;
        A[i][cols - 1] = 50.0;
    }

    for (int j = 0; j < cols; j++) {
        A[0][j]        = 50.0;
        A[rows - 1][j] = 80.0;
    }

    /* ---------- store initial state (k = 0) ---------- */
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            history[idx(0, i, j)] = A[i][j];

    /* ---------- Gauss–Seidel iterations ---------- */
    for (int k = 1; k < nt; k++) {

        for (int i = 1; i < rows - 1; i++) {
            for (int j = 1; j < cols - 1; j++) {
                A[i][j] = 0.25 * (
                    A[i+1][j] + A[i-1][j] +
                    A[i][j+1] + A[i][j-1]
                );
            }
        }

        /* copy FULL mesh into history */
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                history[idx(k, i, j)] = A[i][j];
    }

    /* ---------- write to file ---------- */
    write_history(history);

    return 0;
}

