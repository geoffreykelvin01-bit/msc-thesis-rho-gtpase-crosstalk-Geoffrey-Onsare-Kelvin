# Chapter 4 — ODE Bifurcation Analysis

Scripts for the temporal (non-spatial) Rac–Rho ODE model: phase plane analysis and bifurcation diagrams as the parameter **a5** is varied across the three dynamical scenarios (S1 bistable, S2 oscillatory, S3 coexistence).

## Scripts

### `phase_portrait.m`
Generates multi-panel phase portraits and time series for a chosen scenario (`S1`, `S2`, or `S3`) across three representative `a5` values. For each `a5` value it plots:
- Nullclines (R-nullcline and rho-nullcline) overlaid on a normalised vector field
- Trajectories from one or two initial conditions, integrated with `ode15s`
- Corresponding time series of R(t) and rho(t)

Edit the `scenario` variable at the top of the script to switch between S1/S2/S3. Output figures are saved to a `Nullcline_Analysis_ODE/` folder created relative to the current working directory.

### `plot_matcont_bifurcation.m`
Post-processes MatCont continuation output (`.mat` files exported from MatCont sessions) into publication-style bifurcation diagrams:
- One-parameter continuation: R and rho vs. a5, with stable/unstable branches distinguished and Hopf (H) / limit point (LP) markers labelled
- Two-parameter continuation: Hopf curves in the (a5, a2) plane, with the oscillatory region shaded and Bogdanov–Takens (BT) / Generalized Hopf (GH) points extracted

Requires the corresponding MatCont `.mat` exports (`EP_EP(1).mat`, `EP_EP(2).mat`, `H_H(1)_c1_a2.mat`, etc.) in the working directory.

## Dependencies
MATLAB ODE suite (`ode15s`), MatCont (for generating the `.mat` inputs to `plot_matcont_bifurcation.m`).
