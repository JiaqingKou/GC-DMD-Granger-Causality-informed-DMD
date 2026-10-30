# GC-DMD — VAR(1) benchmark, paper Figures 1–4

Run `main` from MATLAB. The script reproduces the original **separate** figure windows.

- `main.m`: original experiment parameters, data split, GC, lambda sweep, and figure calls.
- `methods/`: DMD, GC-DMD, Oracle and original shared regularized solver (no standalone L1-DMD experiment).
- `figures/plot_fig1.m` ... `plot_fig4.m`: extracted original manuscript plotting code.
- `figures/finalize_paper_figures.m`: **essential original post-plot pass**. Applies mathematical italics, label styling, and endpoint ticks. This step was missing in the previous GitHub packages.
- `utils/`: data generator, diagnostics and original color/style helpers.

**Important:** Do not remove the `finalize_paper_figures(...)` call at the end of `main.m`. The original figure routines deliberately defer some math typography to this final pass.

This package was inspected statically but has not been executed in MATLAB here.

## Unified optimization
`methods/optimize_operator.m` consolidates the original FISTA solver, GC-DMD regularization path and Oracle validation search. Its modes are `solve`, `gc_path`, and `validation`. The original numerical routines and plotting files are otherwise unchanged.

## GC penalty weights
`methods/build_gc_weights.m` replaces `final_gc_weights.m` and
`build_penalty_matrices.m`. It accepts an existing square GC matrix (no
re-estimation) or a snapshot matrix (re-estimates GC for the final refit).
`main.m` reuses its training GC matrix, while the optimizer re-estimates
GC on the combined train/validation segment as before.

## GC intercept specification
The conditional VAR(1) GC regressions include an **explicit intercept** in
**both** the full and reduced models, controlled by `gc_intercept=true` in
`main.m`. For target `i` and candidate driver `j`:

- Full: `x_i(t+1) = b_i + sum_k A_ik x_k(t) + e_F(t)`
- Reduced: `x_i(t+1) = b_i^(R) + sum_(k != j) A_ik^(R) x_k(t) + e_R(t)`

Both are fitted independently by OLS, using **one-step** residuals.
With `T = N_train-1` paired observations and `r` lagged variables, the
full-model residual degrees of freedom are `df2 = T - (r+1)` and the
nested F test has `df1 = 1`. BH-FDR is applied over `r*(r-1)` directed
off-diagonal tests. `info.b_full` exposes the fitted full-model intercepts.

The snapshot demeaning in `main.m` is unchanged to preserve the original
benchmark preprocessing. DMD and the regularized operator optimizer are
**still fitted without an intercept**; only GC regression is changed.
Consequently, results can differ from the prior no-intercept GC version.

The same intercept option is passed into the GC re-estimation on the
combined training+validation segment during the final GC-DMD refit.

**Validation:** static source inspection only; MATLAB runtime not available.


## Automatic EPS figure export

Run `main` in MATLAB. After manuscript formatting, the script exports nine standalone vector EPS files to `results_data/eps/` using `print -depsc2 -painters` (compatible with older MATLAB releases). The original plotting functions are unchanged. Figure 1 consists of three files, Figure 2 of four files, and Figures 3 and 4 of one file each.
