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

model <- fit_cdrgam(
    rt ~ s(word_position, k = 5) +
        irf(surprisal, window = c(0, 2), k = 20),
    impulses = words,
    responses = fixations,
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
    rt ~ irf(surprisal, window = c(0, 2), k = 20),
    impulses = words,
    responses = fixations,
    series = c("subject", "document")
)

model <- fit_cdrgam(design, backend = "sparse")
```

`predict_cdrgam()` accepts new untiled impulse and response streams and rebuilds
their response-level design from the fitted bases and constraints.

## Supported terms

The current implementation supports:

- stationary linear impulse-response functions;
- nonlinear tensor surfaces over lag and impulse value;
- response-aligned `by` modulation;
- lag-by-covariate varying impulse-response functions;
- factor-specific random impulse-response deviations;
- ordinary `mgcv` terms, including response-side random effects; and
- Gaussian and non-Gaussian families through the native `mgcv` backend.

Cubic regression splines are currently the supported IRF basis.

## Fitting backends

`fit_cdrgam()` provides three fitting paths:

- `backend = "mgcv"` fits the compiled response-level design with
  `mgcv::gam()` or `mgcv::bam()`. The result directly inherits from the native
  `gam` or `bam` class.
- `backend = "block"` is an independent dense Gaussian REML reference solver.
- `backend = "sparse"` fits Gaussian identity-link models with sparse penalized
  normal equations and sparse Cholesky factorization. Grouped IRFs remain
  compact until backend assembly.

`backend = "sparse_trust"` selects the sparse engine with its experimental
exact-score, safeguarded trust-region BFGS optimizer. The sparse backends also
support resumable checkpoints and structured progress reporting.

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

## License

`cdrgam` is distributed under the MIT License. See
[cdrgam/LICENSE](cdrgam/LICENSE).
