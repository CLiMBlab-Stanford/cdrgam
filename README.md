# cdrgam

`cdrgam` fits continuous-time deconvolutional regression (CDR) models in R.
It represents impulse-response functions with GAM smooths while compiling
irregular impulse histories into compact response-level designs. This avoids
materializing the much larger history-level design used by a direct `mgcv`
linear-functional fit.

The package is an experimental rewrite. Its API and numerical backends may
change while validation continues.

## Model interface

A model combines ordinary response-aligned GAM terms with `irf()` terms whose
predictors come from an impulse stream:

```r
library(cdrgam)

model <- cdrgam(
    rt ~ s(word_position, k = 5) +
        irf(surprisal, k_l = 20),
    impulses = words,
    responses = fixations,
    window = c(0, 2),
    series = c("subject", "document"),
    impulse_time = "time",
    response_time = "time",
    backend = "sparse"
)
```

Every formula receives a deconvolutional intercept, represented by an implicit
`irf(1)`, unless the formula removes it with `- irf(1)`. This is distinct from
the ordinary response intercept.

Use `prepare_cdrgam()` when the same compiled design will be fitted more than
once:

```r
design <- prepare_cdrgam(
    rt ~ irf(surprisal, k_l = 20),
    impulses = words,
    responses = fixations,
    window = c(0, 2),
    series = c("subject", "document")
)

model <- cdrgam.fit(design, backend = "sparse")
```

`predict()` accepts new untiled streams in a named `newdata` list and rebuilds
their response-level design from the fitted bases and constraints:

```r
prediction <- predict(model, newdata = list(
    impulses = new_words,
    responses = new_fixations
))
```

Set `rescale_predictors = TRUE` to improve numerical conditioning when
continuous predictors or timestamps have awkward units. The transformation
divides by training-data standard deviations but never centers, so interaction
reference points retain their meaning. Binary, categorical, grouping, series,
and constant variables are left alone. History windows are built in native
time units, then lag and time are divided by the same impulse-time standard
deviation. Formulas remain unchanged; new data are transformed with the stored
training statistics, and CDR plots and effect estimates use native source
units.

## Plotting

`plot()` provides CDR-specific views of the fitted terms:

```r
plot(model, view = "irf", select = "surprisal",
     at = list(predictor = c(-1, 0, 1)))
plot(model, view = "predictor", select = "surprisal",
     at = list(lag = c(0.2, 0.5, 1)))
plot(model, view = "surface", select = "surprisal")
```

Grouped terms support population, deviation, and conditional views. CDR plots
show additive effects on the linear-predictor scale and can include conditional
or smoothing-parameter-adjusted uncertainty. `save_cdrgam_plots()` writes the
same panels to PDF or raster devices. Native `mgcv` fits retain access to the
translated GAM display through `plot(model, view = "gam")` or
`mgcv::plot.gam(as_gam(model))`. The first form restores scaled ordinary-smooth
axes to source units; the explicit `as_gam()` escape hatch exposes the literal
internal mgcv parameterization.

## Supported terms

The current implementation supports:

- stationary linear and nonlinear impulse-response functions;
- multi-predictor tensor interactions with mixed linear and smooth axes;
- response-aligned `by` modulation;
- nonstationary response-time axes;
- grouped impulse-response deviations;
- ordinary `mgcv` terms, including response-side random effects; and
- Gaussian and non-Gaussian families through the native `mgcv` backend.

`k_l`, `k_t`, and `k_p` control the lag, response-time, and predictor basis
dimensions. `NULL` predictor entries are linear; use an R list such as
`k_p = list(NULL, 4)` when an interaction mixes linear and smooth predictors.
The corresponding `bs_l`, `bs_t`, and `bs_p` arguments accept standard numeric
`mgcv` marginal bases. Smooth non-lag tensor marginals are centered to separate
an interaction from its lower-order effects.

## Fitting backends

`cdrgam()` provides three fitting paths:

- `backend = "mgcv"` fits the compiled response-level design with
  `mgcv::gam()` or `mgcv::bam()`. The result directly inherits from the native
  `gam` or `bam` class.
- `backend = "block"` is an independent dense Gaussian REML reference solver.
- `backend = "sparse"` fits Gaussian identity-link models with sparse penalized
  normal equations and sparse Cholesky factorization. Grouped IRFs remain
  compact until backend assembly.

The sparse backend's experimental exact-score, safeguarded trust-region BFGS
optimizer is selected with
`sparse_control = list(gradient = "exact", outer_optimizer = "bfgs_trust")`.
Its checkpoints preserve
the complete optimizer state, including BFGS curvature and trust radius, so an
interrupted run continues its original trajectory. If the trust radius reaches
the numerical floor during fitting, the resolved Hessian strategy is used to
certify practical convergence in identifiable curvature directions or reset
curvature for a bounded recovery attempt. The default automatic strategy uses
the analytic Hessian when it fits the process memory budget and otherwise
differences exact gradients. Near-neutral steps are accepted only when their
exact score reaches tolerance or improves, which prevents factorization-level
objective noise from collapsing the trust radius.

With `sparse_control = list(schur = "always")`, response-aligned random-effect
blocks are eliminated through a consolidated Schur transfer matrix. The
factorization, exact-score trace, and multi-right-hand-side solves use batched
BLAS operations. A threaded BLAS implementation such as OpenBLAS materially
reduces runtime for large Schur cores; R's reference BLAS remains supported
but executes these kernels serially.

Set `sparse_control$boundary_action = "reduce"` to convert certified IRF
curvature boundaries into explicit reduced bases and reoptimize the remaining
smoothing parameters. The user formula is preserved, and the fitted object
records the removed penalty subspaces, original and effective dimensions, and
the conditional status of subsequent inference. The default `"report"`
records the same diagnosis without altering the fitted representation. Reports
name the affected component (for example, lag curvature, predictor curvature,
or group-deviation magnitude). A full-rank boundary would erase an entire
term, so it is reported for confirmation and is never removed automatically.

The independent backends are used to test compilation, fitting, prediction,
and inference against one another.

## Installation

Install from the repository root:

```sh
R CMD INSTALL cdrgam
```

The package requires R, `Matrix`, `methods`, and `mgcv`. It compiles a small C
extension for the optional Schur-complement factorization path.

## Tests

Run the standalone test scripts from the repository root:

```sh
for test_file in cdrgam/tests/*.R; do
    Rscript "$test_file" || exit 1
done
```

The suite covers compressed-design equivalence, formula compilation,
identifiability, nonlinear and varying effects, mixed models, sparse
factorization, prediction, checkpoint recovery, optimizer behavior, and
synthetic recovery.

Longer empirical and performance checks live under `validation/`. Their
generated outputs and benchmark reports are intentionally excluded from Git.

## Development

Development takes place on `dev` or a feature branch. See
[CONTRIBUTING.md](CONTRIBUTING.md) for testing, commit attribution, and release
requirements. Instructions for AI coding agents are in [AGENTS.md](AGENTS.md),
and repository-level AI assistance is recorded in
[AI_PROVENANCE.md](AI_PROVENANCE.md).

## AI-assisted development

This codebase was developed almost entirely with AI coding assistance. Human
maintainers selected the requirements, reviewed the generated code, ran the
applicable tests, and accept responsibility for the published result.

The tools and models used are recorded in
[AI_PROVENANCE.md](AI_PROVENANCE.md). Where available, individual commits also
contain `Assisted-by:` trailers.

## License

`cdrgam` is distributed under the MIT License. See
[cdrgam/LICENSE](cdrgam/LICENSE).
