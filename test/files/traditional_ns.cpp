/*
 * traditional_ns.cpp — Traditional Monolithic Navier-Stokes Solver
 * =================================================================
 *
 * 2D Lid-Driven Cavity, Chorin's Projection Method
 * No tiling, no game engine techniques. Pure monolithic approach.
 *
 * Provides two variants:
 *   - Jacobi pressure solver (100 iterations)
 *   - Multigrid pressure solver (V-cycle)
 *
 * Compile (standalone test):
 *   g++ -O3 -fopenmp -march=native -DTRADITIONAL_STANDALONE \
 *       -o traditional_ns traditional_ns.cpp -lm
 *
 * Used by benchmark.cpp when compiled without -DTRADITIONAL_STANDALONE
 */

#ifndef TRADITIONAL_NS_H
#define TRADITIONAL_NS_H

#include "ns_common.h"

// =============================================================================
// Full NS step (monolithic — operates on entire grid)
// =============================================================================
static void ns_step_monolithic(Field& u, Field& v, Field& p,
                               Field& u_lap, Field& v_lap,
                               Field& dudx, Field& dudy,
                               Field& dvdx, Field& dvdy,
                               Field& u_star, Field& v_star,
                               Field& div_f, Field& rhs,
                               double dx, double dt, double nu,
                               bool use_mg)
{
    int ny = u.ny, nx = u.nx;

    // Step 1: Advection + Diffusion → u*
    compute_laplacian(u, u_lap, dx);
    compute_laplacian(v, v_lap, dx);
    compute_ddx(u, dudx, dx); compute_ddy(u, dudy, dx);
    compute_ddx(v, dvdx, dx); compute_ddy(v, dvdy, dx);

    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            double adv_u = u(i,j)*dudx(i,j) + v(i,j)*dudy(i,j);
            double adv_v = u(i,j)*dvdx(i,j) + v(i,j)*dvdy(i,j);
            u_star(i,j) = u(i,j) + dt * (-adv_u + nu * u_lap(i,j));
            v_star(i,j) = v(i,j) + dt * (-adv_v + nu * v_lap(i,j));
        }
    }
    apply_velocity_bc(u_star, v_star);

    // Step 2: Pressure Poisson (∇²p = (1/dt)∇·u*)
    compute_ddx(u_star, dudx, dx);
    compute_ddy(v_star, dvdy, dx);
    double inv_dt = 1.0 / dt;
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++)
            rhs(i,j) = (dudx(i,j) + dvdy(i,j)) * inv_dt;

    if (use_mg)
        pressure_multigrid(p, rhs, dx);
    else
        pressure_jacobi(p, rhs, dx, JACOBI_ITERS);

    // Step 3: Projection (u = u* - dt·∇p)
    compute_ddx(p, dudx, dx);
    compute_ddy(p, dvdy, dx);
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            u(i,j) = u_star(i,j) - dt * dudx(i,j);
            v(i,j) = v_star(i,j) - dt * dvdy(i,j);
        }
    }
    apply_velocity_bc(u, v);
}

// =============================================================================
// Run monolithic solver for N timesteps
// =============================================================================
struct TraditionalResult {
    double time;
    // Fields are modified in-place
};

static TraditionalResult run_monolithic(Field& u, Field& v, Field& p,
                                        int n_steps, bool use_mg)
{
    int n = GRID_SIZE;
    Field u_lap(n,n), v_lap(n,n), dudx(n,n), dudy(n,n);
    Field dvdx(n,n), dvdy(n,n), u_star(n,n), v_star(n,n);
    Field div_f(n,n), rhs(n,n);

    Timer tm; tm.start();
    for (int s = 0; s < n_steps; s++) {
        ns_step_monolithic(u, v, p, u_lap, v_lap, dudx, dudy, dvdx, dvdy,
                           u_star, v_star, div_f, rhs, DX, DT, NU, use_mg);
        if ((s+1) % 20 == 0)
            printf("    step %d/%d\n", s+1, n_steps);
    }
    double elapsed = tm.elapsed();

    u_lap.free_mem(); v_lap.free_mem(); dudx.free_mem(); dudy.free_mem();
    dvdx.free_mem(); dvdy.free_mem(); u_star.free_mem(); v_star.free_mem();
    div_f.free_mem(); rhs.free_mem();

    return {elapsed};
}

// =============================================================================
// Standalone test
// =============================================================================
#ifdef TRADITIONAL_STANDALONE
int main() {
    int nt = omp_get_max_threads();
    printf("Traditional NS Solver (Monolithic)\n");
    printf("Grid: %d×%d, Steps: %d, Threads: %d\n\n", GRID_SIZE, GRID_SIZE, NUM_TIMESTEPS, nt);

    Field u, v, p;
    init_fields(u, v, p);

    printf("[1] Jacobi (%d iterations):\n", JACOBI_ITERS);
    Field u1(GRID_SIZE,GRID_SIZE), v1(GRID_SIZE,GRID_SIZE), p1(GRID_SIZE,GRID_SIZE);
    u1.copy_from(u); v1.copy_from(v); p1.copy_from(p);
    auto r1 = run_monolithic(u1, v1, p1, NUM_TIMESTEPS, false);
    printf("    -> %.2f s\n\n", r1.time);

    printf("[2] Multigrid (V%d, pre%d/post%d):\n", MG_VCYCLES, MG_PRE_SMOOTH, MG_POST_SMOOTH);
    Field u2(GRID_SIZE,GRID_SIZE), v2(GRID_SIZE,GRID_SIZE), p2(GRID_SIZE,GRID_SIZE);
    u2.copy_from(u); v2.copy_from(v); p2.copy_from(p);
    auto r2 = run_monolithic(u2, v2, p2, NUM_TIMESTEPS, true);
    printf("    -> %.2f s (%.2fx faster)\n", r2.time, r1.time / r2.time);

    u.free_mem(); v.free_mem(); p.free_mem();
    u1.free_mem(); v1.free_mem(); p1.free_mem();
    u2.free_mem(); v2.free_mem(); p2.free_mem();
    return 0;
}
#endif

#endif // TRADITIONAL_NS_H
