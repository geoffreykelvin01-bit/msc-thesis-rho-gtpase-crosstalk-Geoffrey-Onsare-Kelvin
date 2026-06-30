# Chapter 5 — Diffusion Analysis (Local Perturbation Analysis)

Scripts for Local Perturbation Analysis (LPA) of the spatial Rac–Rho model: dispersion relations, averaged Jacobian stability analysis, and combined LPA/ODE bifurcation visualisation.

## Scripts

### `Averaged_Jacobian.m`
Computes the averaged (spatially homogenised) Jacobian for the LPA system and performs dispersion relation / stability analysis across the diffusion coefficient. Identifies zero-crossings to flag critical diffusion values, shades unstable wavenumber bands, and exports per-scenario (S1/S2/S3) figures into organised subfolders.

### `combined_lpa_ode_plotter.m`
Comprehensive plotting suite combining one- and two-parameter bifurcation diagrams for both the ODE and LPA systems side by side. Includes helper functions for loading and indexing MatCont data matrices across file types (`EP_EP`, `H_LC`, `LP_LP`, `H_H`, `BP_BB`) for both ODE and LPA continuations.

Expects MatCont `.mat` exports under an `OSCILLATORY_trial/` folder (created relative to the current working directory), with `LPA<suffix>/` and `ODE<suffix>/` subfolders for the LPA and ODE continuation data respectively. Output figures are written to `OSCILLATORY_trial/figures_for_<set_num>/`.

## Dependencies
MatCont (for `.mat` continuation data), MATLAB base plotting/linear algebra functions.

## Notes
Row index constants for MatCont data matrices differ between ODE and LPA layouts and between one- and two-parameter continuations — see in-script comments/constants before adapting to new data.
