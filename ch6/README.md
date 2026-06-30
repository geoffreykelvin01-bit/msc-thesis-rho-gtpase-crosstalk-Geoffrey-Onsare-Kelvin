# Chapter 6 — Numerical Simulations (PDE)

Full reaction-diffusion PDE simulations of the four-component Rac–Rho model (active/inactive Rac, active/inactive Rho), solved with MATLAB's `pdepe`.

## Scripts

### `PDEPE_SOLVER.m`
End-to-end PDE simulation and figure-generation pipeline for a chosen dynamical regime (`b` bistable, `o` oscillatory, `c` coexistence):
- Solves the system to steady state (via `fsolve`/`ode45` fallback), then perturbs it with a configurable initial condition (`manual`, `random`, `cosine`, or `local` bump)
- Integrates the PDE with `pdepe` for each initial condition in `IC_list`
- Exports spatial profile plots, space–time surface (`pcolor`) plots, L² norm convergence plots, mass-conservation checks, initial-condition profiles, and multi-panel combined figures
- All figures exported as high-DPI PNGs into a regime/perturbation-labelled subfolder

Edit `regime`, `pert_type`, `k_mode`, `pert_mag`, `t_end`, and `IC_list` at the top of the script to configure a run.

## Dependencies
MATLAB PDE toolbox functions (`pdepe`), Optimization Toolbox (`fsolve`).
