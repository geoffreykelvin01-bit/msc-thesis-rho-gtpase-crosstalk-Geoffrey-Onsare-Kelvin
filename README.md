# Rho GTPase Crosstalk: Reaction-Diffusion Models of Cell Polarisation
**Masters Thesis Computational Code — Geoffrey Onsare Kelvin**
Department of Mathematics | Mathematical Biology

---

## Overview

This repository contains all MATLAB code developed for a Masters thesis investigating the spatiotemporal dynamics of Rho GTPase crosstalk (Rac and Rho) in cell migration and polarisation. The work combines ODE bifurcation analysis, Local Perturbation Analysis (LPA), PDE numerical simulations, and dispersion relation analysis to study how mutual inhibition between Rac and Rho gives rise to cell polarity.

The modelling framework is grounded in experimental data from Nanda et al. (2023), with the dimensionless parameter **a5** (Rac-to-Rho feedforward strength, mediated by Arhgef11/12) serving as the primary bifurcation parameter across three dynamical scenarios:

| Scenario | Description |
|----------|-------------|
| **S1** | Bistable — two stable steady states, wave-pinning polarisation possible |
| **S2** | Oscillatory — limit cycle behaviour |
| **S3** | Coexistence — bistable and oscillatory regions overlap |

---

## Repository Structure

```
.
├── ch4/    # ODE bifurcation analysis (phase portraits, MatCont bifurcation diagrams)
├── ch5/    # Diffusion analysis / Local Perturbation Analysis (LPA, dispersion relations)
├── ch6/    # PDE numerical simulations (full reaction-diffusion solver, pdepe)
└── README.md
```

Each chapter folder has its own `README.md` describing the scripts inside it, their inputs/outputs, and dependencies.

### ch4 — ODE Bifurcation Analysis
- `phase_portrait.m` — nullcline/phase portraits and time series across S1/S2/S3
- `plot_matcont_bifurcation.m` — one- and two-parameter bifurcation diagrams from MatCont output

### ch5 — Diffusion Analysis (LPA)
- `Averaged_Jacobian.m` — averaged Jacobian / dispersion relation stability analysis
- `combined_lpa_ode_plotter.m` — combined LPA + ODE bifurcation plotting suite

### ch6 — Numerical Simulations (PDE)
- `PDEPE_SOLVER.m` — full four-component reaction-diffusion PDE solver and figure pipeline

---

## Requirements

- MATLAB (developed/tested on recent releases)
- [MatCont](https://sourceforge.net/projects/matcont/) for continuation analysis (ch4, ch5)
- Optimization Toolbox (`fsolve`) for ch6

## Citation

If you use this code, please cite the associated Masters thesis (details to be added upon submission).
