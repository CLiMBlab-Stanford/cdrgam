# Repository instructions for AI agents

## Required reading

Before making a material change, read `AI_PROVENANCE.md` and register
the provider, tool or agent, exposed model identifier, period, and roles
under **Systems used**. Update an existing row when appropriate instead
of adding one for every session.

Before creating or revising documentation, docstrings, or explanatory
code comments, read and follow `WRITING_POLICY.md`. Apply it to prose
touched by the task without rewriting unrelated text.

Use only identity information exposed by the environment. Never infer a
model version or store prompts, secrets, private context, conversation
logs, or session identifiers.

## Repository structure

The R package occupies the repository root:

- `R/formula.R` parses
  [`irf()`](https://climblab-stanford.github.io/cdrgam/dev/reference/irf.md)
  terms and compiles stream histories.
- `R/compressed.R` integrates compact designs with native `mgcv` fits.
- `R/block.R` implements the dense Gaussian REML reference backend.
- `R/sparse.R` implements the scalable sparse Gaussian REML backend.
- `R/solver-control.R` contains checkpoints, progress reporting,
  numerical Hessians, and the experimental trust optimizer.
- `R/prediction.R` rebuilds designs for new streams.
- `src/` contains registered native routines.
- `tests/` contains standalone regression and recovery tests.
- `validation/` contains longer empirical and performance experiments.

The historical implementation under `archive/` is reference-only,
ignored by Git, and not authoritative for the active API or model
semantics.

## Correctness and validation

Check mathematical and API claims against the implementation. Preserve
the independent backend comparisons: native `mgcv`, dense block REML,
and sparse REML provide complementary checks on preprocessing,
optimization, prediction, and inference.

Run the narrowest relevant tests while developing, then run all
standalone tests before handing off a broad change:

``` sh
for test_file in tests/*.R; do
    Rscript "$test_file" || exit 1
done
```

Use a clean temporary library when installed-package behavior or
namespace contents matter. Do not commit installed libraries, compiled
objects, package tarballs, checkpoints, logs, plots, fitted objects, or
generated validation results.

Long-running validation may have active cluster jobs and resumable
checkpoints. Inspect job and checkpoint state before changing optimizer
or serialization contracts. Do not cancel, overwrite, or invalidate a
run unless the user asks.

## Commits and publication

Commit, push, tag, or create or update a pull request only when the user
explicitly requests that action. Development belongs on `dev` or a
feature branch. Changes enter the release-only `main` branch through a
pull request after its hosted release gate and cross-platform checks
pass.

Before committing, verify that the existing Git configuration provides
both the requesting human’s `user.name` and `user.email`. Before using a
hosting CLI or API, verify that the authenticated session belongs to
that human. Do not change identity settings, use `--author`, or silently
substitute an agent, bot, service account, or another person.

If identity is missing or ambiguous, stop before publication. Leave the
work ready and provide the proposed commit message and applicable
`Assisted-by:` trailers. Never request, print, store, or copy
authentication secrets into the repository.

AI systems must not be Git authors, copyright holders, DCO signatories,
or `Signed-off-by:` identities. The human contributor owns and approves
the change. Follow the attribution and signing rules in
`CONTRIBUTING.md`.

Do not fabricate a cryptographic signature. Use a configured human
signing mechanism only with authorization. Add `Signed-off-by:` only
when the human has authorized the certification it represents.

## Writing

Write for researchers and engineers who may know GAMs but not this
implementation. State observable behavior before implementation detail.
Use the established terms *impulse stream*, *response stream*, *history
link*, *impulse-response function*, *response-level design*, *smoothing
parameter*, and *backend* consistently.

Document constraints and experimental status accurately. In particular,
do not describe custom Gaussian backends as native `mgcv` objects, do
not generalize Gaussian-only features to every family, and do not claim
convergence merely because an optimizer returned a finite fit.

AI disclosure belongs in `AI_PROVENANCE.md` and commit or pull-request
metadata, not in ordinary technical documentation.
