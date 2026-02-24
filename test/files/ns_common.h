/*
 * ns_common.h — Shared types, parameters, and core numerics
 * ==========================================================
 *
 * Used by both traditional_ns.cpp and game_engine_ns.cpp
 */

#ifndef NS_COMMON_H
#define NS_COMMON_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>
#include <chrono>
#include <omp.h>

// =============================================================================
// Parameters
// =============================================================================
static const int    GRID_SIZE         = 2048;
static const int    NUM_TIMESTEPS     = 60;
static const int    TILE_SIZE         = 128;
static const int    HALO              = 1;
static const int    TILE_FULL         = TILE_SIZE + 2 * HALO;  // 130

static const double RE                = 100.0;
static const double LID_VELOCITY      = 1.0;
static const double DX                = 1.0 / (GRID_SIZE - 1);
static const double NU                = LID_VELOCITY / RE;
static const double DT                = 5e-5;

// Pressure solver
static const int    JACOBI_ITERS      = 100;

// Multigrid
static const int    MG_LEVELS         = 3;
static const int    MG_PRE_SMOOTH     = 3;
static const int    MG_POST_SMOOTH    = 3;
static const int    MG_VCYCLES        = 2;

// Culling / reuse
static const double CULL_THRESHOLD    = 1e-6;
static const double REUSE_THRESHOLD   = 1e-5;

// Tiles per direction
static const int    NTILES_DIR        = GRID_SIZE / TILE_SIZE;
static const int    NTILES            = NTILES_DIR * NTILES_DIR;

// =============================================================================
// Timer
// =============================================================================
struct Timer {
    std::chrono::high_resolution_clock::time_point t0;
    void start() { t0 = std::chrono::high_resolution_clock::now(); }
    double elapsed() const {
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(t1 - t0).count();
    }
};

// =============================================================================
// 2D array wrapper (contiguous, row-major)
// =============================================================================
struct Field {
    double* data;
    int ny, nx;

    Field() : data(nullptr), ny(0), nx(0) {}
    Field(int ny_, int nx_) : ny(ny_), nx(nx_) {
        data = (double*)aligned_alloc(64, ny * nx * sizeof(double));
        memset(data, 0, ny * nx * sizeof(double));
    }
    ~Field() { /* Manual free — we manage lifetime explicitly */ }

    void alloc(int ny_, int nx_) {
        ny = ny_; nx = nx_;
        data = (double*)aligned_alloc(64, ny * nx * sizeof(double));
        memset(data, 0, ny * nx * sizeof(double));
    }
    void free_mem() { if (data) { free(data); data = nullptr; } }
    void zero() { memset(data, 0, ny * nx * sizeof(double)); }
    void copy_from(const Field& src) {
        memcpy(data, src.data, ny * nx * sizeof(double));
    }

    inline double& operator()(int i, int j) { return data[i * nx + j]; }
    inline double  operator()(int i, int j) const { return data[i * nx + j]; }

    // Periodic access
    inline double at(int i, int j) const {
        int ii = (i + ny) % ny;
        int jj = (j + nx) % nx;
        return data[ii * nx + jj];
    }
};

// =============================================================================
// Boundary conditions
// =============================================================================
static inline void apply_velocity_bc(Field& u, Field& v) {
    int ny = u.ny, nx = u.nx;
    for (int j = 0; j < nx; j++) {
        u(0, j)    = 0.0;
        u(ny-1, j) = LID_VELOCITY;
        v(0, j)    = 0.0;
        v(ny-1, j) = 0.0;
    }
    for (int i = 0; i < ny; i++) {
        u(i, 0)    = 0.0;
        u(i, nx-1) = 0.0;
        v(i, 0)    = 0.0;
        v(i, nx-1) = 0.0;
    }
}

static inline void apply_pressure_bc(Field& p) {
    int ny = p.ny, nx = p.nx;
    for (int j = 0; j < nx; j++) {
        p(0, j)    = p(1, j);
        p(ny-1, j) = p(ny-2, j);
    }
    for (int i = 0; i < ny; i++) {
        p(i, 0)    = p(i, 1);
        p(i, nx-1) = p(i, nx-2);
    }
}

// =============================================================================
// Core differential operators
// =============================================================================
static inline void compute_laplacian(const Field& f, Field& lap, double dx) {
    double inv_dx2 = 1.0 / (dx * dx);
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < f.ny; i++)
        for (int j = 0; j < f.nx; j++)
            lap(i, j) = (f.at(i-1, j) + f.at(i+1, j) +
                         f.at(i, j-1) + f.at(i, j+1) - 4.0 * f(i, j)) * inv_dx2;
}

static inline void compute_ddx(const Field& f, Field& out, double dx) {
    double inv_2dx = 1.0 / (2.0 * dx);
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < f.ny; i++)
        for (int j = 0; j < f.nx; j++)
            out(i, j) = (f.at(i, j+1) - f.at(i, j-1)) * inv_2dx;
}

static inline void compute_ddy(const Field& f, Field& out, double dx) {
    double inv_2dx = 1.0 / (2.0 * dx);
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < f.ny; i++)
        for (int j = 0; j < f.nx; j++)
            out(i, j) = (f.at(i+1, j) - f.at(i-1, j)) * inv_2dx;
}

// =============================================================================
// Jacobi pressure solver
// =============================================================================
static inline void pressure_jacobi(Field& p, const Field& rhs, double dx, int n_iters) {
    double dx2 = dx * dx;
    Field p_new(p.ny, p.nx);
    for (int iter = 0; iter < n_iters; iter++) {
        #pragma omp parallel for schedule(static) collapse(2)
        for (int i = 0; i < p.ny; i++)
            for (int j = 0; j < p.nx; j++)
                p_new(i, j) = 0.25 * (p.at(i-1, j) + p.at(i+1, j) +
                                       p.at(i, j-1) + p.at(i, j+1) -
                                       dx2 * rhs(i, j));
        apply_pressure_bc(p_new);
        double* tmp = p.data; p.data = p_new.data; p_new.data = tmp;
    }
    p_new.free_mem();
}

// =============================================================================
// Multigrid V-cycle pressure solver
// =============================================================================
static inline void mg_restrict(const Field& fine, Field& coarse) {
    int cny = coarse.ny, cnx = coarse.nx;
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < cny; i++)
        for (int j = 0; j < cnx; j++) {
            int fi = 2 * i, fj = 2 * j;
            if (fi + 1 < fine.ny && fj + 1 < fine.nx)
                coarse(i, j) = 0.25 * (fine(fi, fj) + fine(fi+1, fj) +
                                        fine(fi, fj+1) + fine(fi+1, fj+1));
            else
                coarse(i, j) = fine(std::min(fi, fine.ny-1),
                                    std::min(fj, fine.nx-1));
        }
}

static inline void mg_prolongate(const Field& coarse, Field& fine) {
    int cny = coarse.ny, cnx = coarse.nx;
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < fine.ny; i++)
        for (int j = 0; j < fine.nx; j++)
            fine(i, j) = coarse(std::min(i / 2, cny - 1),
                                std::min(j / 2, cnx - 1));
}

static inline void mg_residual(const Field& p, const Field& rhs, Field& r, double dx) {
    double inv_dx2 = 1.0 / (dx * dx);
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < p.ny; i++)
        for (int j = 0; j < p.nx; j++) {
            double Lp = (p.at(i-1, j) + p.at(i+1, j) +
                         p.at(i, j-1) + p.at(i, j+1) - 4.0 * p(i, j)) * inv_dx2;
            r(i, j) = rhs(i, j) - Lp;
        }
}

static inline void mg_smooth(Field& p, const Field& rhs, double dx, int n_iters) {
    double dx2 = dx * dx;
    Field tmp(p.ny, p.nx);
    for (int iter = 0; iter < n_iters; iter++) {
        #pragma omp parallel for schedule(static) collapse(2)
        for (int i = 0; i < p.ny; i++)
            for (int j = 0; j < p.nx; j++)
                tmp(i, j) = 0.25 * (p.at(i-1, j) + p.at(i+1, j) +
                                     p.at(i, j-1) + p.at(i, j+1) -
                                     dx2 * rhs(i, j));
        apply_pressure_bc(tmp);
        double* sw = p.data; p.data = tmp.data; tmp.data = sw;
    }
    tmp.free_mem();
}

static void mg_vcycle(Field& p, const Field& rhs, double dx,
                      int level, int max_level)
{
    mg_smooth(p, rhs, dx, MG_PRE_SMOOTH);
    if (level < max_level && std::min(p.ny, p.nx) >= 8) {
        Field r(p.ny, p.nx);
        mg_residual(p, rhs, r, dx);
        int cny = p.ny / 2, cnx = p.nx / 2;
        Field r_c(cny, cnx);
        Field e_c(cny, cnx);
        mg_restrict(r, r_c);
        mg_vcycle(e_c, r_c, dx * 2.0, level + 1, max_level);
        Field correction(p.ny, p.nx);
        mg_prolongate(e_c, correction);
        #pragma omp parallel for schedule(static) collapse(2)
        for (int i = 0; i < p.ny; i++)
            for (int j = 0; j < p.nx; j++)
                p(i, j) += correction(i, j);
        r.free_mem(); r_c.free_mem(); e_c.free_mem(); correction.free_mem();
    }
    mg_smooth(p, rhs, dx, MG_POST_SMOOTH);
}

static inline void pressure_multigrid(Field& p, const Field& rhs, double dx) {
    for (int v = 0; v < MG_VCYCLES; v++)
        mg_vcycle(p, rhs, dx, 0, MG_LEVELS);
}

// =============================================================================
// Init fields
// =============================================================================
static inline void init_fields(Field& u, Field& v, Field& p) {
    u.alloc(GRID_SIZE, GRID_SIZE);
    v.alloc(GRID_SIZE, GRID_SIZE);
    p.alloc(GRID_SIZE, GRID_SIZE);
    for (int j = 0; j < GRID_SIZE; j++)
        u(GRID_SIZE - 1, j) = LID_VELOCITY;
}

#endif // NS_COMMON_H
