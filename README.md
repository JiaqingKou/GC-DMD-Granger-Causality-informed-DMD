# GC-DMD — MATLAB manuscript reproduction

# Please cite this reference:
# Yihua Pan, & Jiaqing Kou. (2026). Causality-informed dynamic mode decomposition for unsteady flow analysis. arXiv preprint arXiv:2610.39522.

## Run and figure mapping

|Entry|Case|Manuscript figures|Outputs|
|-|-|-|-|
|`main.m`|VAR(1)|Figures 1–4|EPS|
|`main\_2.m`|Cavity SVD/POD comparison|Figure 5|EPS / PNG|
|`main\_3.m`|Cavity Re=13000|Figures 6–7|EPS / PNG|
|`main\_4.m`|Cavity Re=19000|Figures 8–12|EPS / PNG|
|`main\_5.m`|Cavity Re=30000|Appendix B1–B2|EPS / PNG|
|`main\_6.m`|Mach 0.9 jet|Figures 13–16|EPS / PNG|
|`main\_7.m`|L1-DMD benchmark|Appendix A1–A2|EPS|

Run **one** `main\*.m` script at a time from MATLAB. Every script initializes its own paths and saves only its manuscript figures to `results/`. Re-running one main replaces only that main's figures, independently of execution order.

## Directory structure

* `methods/`: numerical methods (snapshot decomposition, DMD, conditional GC, FISTA/weighted L1, modal analysis).
* `utils/`: dataset loading, case-level diagnostics, plotting styles, and output-directory management.
* `figures/`: original manuscript figure drawing functions, grouped by manuscript figure number; a function may create several MATLAB windows.
* `data/`: input MAT files; `operator\_bank\_20.mat` is included, large CFD datasets must be supplied separately.
* `results/`: generated EPS and PNG images only.

## Main script reading order

`0. Settings → 1. Load snapshots → 2. Snapshot SVD/POD (when applicable) → 3. DMD → 4. GC → 5. GC-DMD → 6+. Modal analysis, forecasts, sparsity, and manuscript figures.`

The important algorithm calls are **visible in the main scripts** rather than hidden behind an all-in-one identification wrapper. The cavity and jet cases use an explicit chronological training/validation split; final DMD and GC-DMD operators are refitted on the modeling interval. GC full/reduced regressions include an intercept; BH-FDR q=0.05; normalized GC rowsum; eta=1; diagonal coefficients are not penalized.

## Numerical fidelity

The supplied project uses `compute\_svd.m` as the naming convention for its **snapshot-correlation eigendecomposition**. This name does not imply calling MATLAB's `svd()` directly. The cavity scripts retain the original QR re-orthogonalization and truncation rules. The existing `methods/optimize\_operator.m` and the figure drawing function bodies are retained; the refactor changes only the *visibility and organization* of the identification sequence in the main scripts.

The original 20-trial Appendix A experiment uses one ground-truth operator with 20 noise seeds, producing 20 separately identified operators. The code does not claim 20 distinct ground-truth systems.

## Data files

Place `Cavity13k.mat`, `Cavity19k.mat`, `Cavity30k.mat`, and `jetLES.mat` under `data/` as needed. Run a case only after its dataset is available. The package includes `operator\_bank\_20.mat` for Appendix A.

## Validation status

Static source checks only; full MATLAB numerical and figure reproduction requires the datasets and MATLAB runtime. No image or figure body was intentionally redesigned.

