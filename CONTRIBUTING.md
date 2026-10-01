# Contributing

## Development workflow

Develop changes on `dev` or a focused feature branch. Keep `main` for reviewed
releases. Make commits small enough to review and include tests for changed
behavior.

Install the package and run the standalone tests from the repository root:

```sh
R CMD INSTALL .

for test_file in tests/*.R; do
    Rscript "$test_file" || exit 1
done
```

During development, use `scripts/test` to install the current source into a
temporary library and run one focused scope:

```sh
./scripts/test api
./scripts/test design
./scripts/test gaussian
./scripts/test generalized
./scripts/test distributional
./scripts/test optimizer
./scripts/test recovery
```

The script also accepts the stem of one file in `tests/` or `all`.

Use a clean temporary R library when testing package installation or namespace
behavior. Run relevant scripts under `validation/` for changes to numerical
methods, scaling behavior, recovery, or empirical model results. Generated
validation results are not committed.

Downstream packages should use exported functions rather than fitted-object
implementation fields. Extend `fit_metadata()`, `fit_report()`, or
`predict_components()` when a downstream tool needs additional stable model
information. Test changes to these functions against native `mgcv`, block, and
sparse fits when the reported field differs by backend.

Before proposing a release, update `Version` in `DESCRIPTION` and document
intentional incompatibilities and migration steps. The hosted release gate is
the authoritative clean-environment check.

## AI-assisted contributions

Disclose material AI assistance in every affected commit by adding one trailer
for each tool and model combination:

```text
Assisted-by: <tool-or-agent>:<model-identifier>
```

Use the most specific model identifier exposed by the tool. If the underlying
model is not disclosed, use the visible product or tool label without guessing.
Place trailers after a blank line at the end of the commit message.

The human contributor remains the commit author and is responsible for review
and testing. Do not list an AI system in `Signed-off-by:` and do not rewrite
published history solely to add missing trailers.

## Commit and pull-request identity

An agent may create a commit or pull request only when the requesting human has
explicitly authorized it. A commit may use the human's name and email only when
both are already available from Git configuration. A pull request may use the
human's account only when the authenticated hosting session clearly belongs to
that person.

An agent must not:

- invent or set a name or email to make a commit appear human-authored;
- use `--author` to override missing or ambiguous identity;
- silently substitute an agent, bot, service account, or another person;
- conceal bot or service attribution when its use was explicitly authorized;
  or
- request, print, store, or copy authentication secrets into the repository.

If identity or authentication is unavailable or ambiguous, stop before the
publication action. Leave the changes ready and provide the proposed commit or
pull-request metadata.

Cryptographic signing must use the human's configured signing mechanism. An
agent must never fabricate a signature. Add `Signed-off-by:` only when the
human has authorized the certification it represents.

## Documentation

Documentation, docstrings, and explanatory comments must follow
`WRITING_POLICY.md`. Check technical claims against the implementation and use
the repository's established terminology.

Keep AI disclosure in `AI_PROVENANCE.md`, commit trailers, and pull-request
metadata instead of ordinary technical documentation.

## Releases and compatibility

`main` contains released code. Prepare releases on `dev` or a release branch
and merge them into `main` only through a pull request. Every pull request to
`main` must change `Version` in `DESCRIPTION` to a later
`MAJOR.MINOR.PATCH` value. Repository protection requires the release gate and
the Linux, macOS, and Windows package checks to pass before merge. The release
gate runs the complete test suite through `R CMD check` and builds the complete
pkgdown site.

Use patch releases for compatible fixes. During the 0.x series, use minor
releases for new features and intentional interface changes. Retain older
behavior when the benefit is clear and the implementation remains readable and
inexpensive to maintain or run.

After the pull request merges, create an annotated `vMAJOR.MINOR.PATCH` tag on
the validated merge commit, using a concise human-written release summary as
the tag message. The tag workflow rejects versions that do not match
`DESCRIPTION`, lightweight or empty tags, and commits outside `main`, then
uses that message as the corresponding GitHub Release description. Verify the
tag workflow and final description before treating publication as complete.
Never move or replace a published tag.
