# Contributing

## Development workflow

Develop changes on `dev` or a focused feature branch. Keep `main` for reviewed
releases. Make commits small enough to review and include tests for changed
behavior.

Install the package and run the standalone tests from the repository root:

```sh
R CMD INSTALL cdrgam

for test_file in cdrgam/tests/*.R; do
    Rscript "$test_file" || exit 1
done
```

Use a clean temporary R library when testing package installation or namespace
behavior. Run relevant scripts under `validation/` for changes to numerical
methods, scaling behavior, recovery, or empirical model results. Generated
validation results are not committed.

Before merging a release, update `Version` in `cdrgam/DESCRIPTION`, run the
complete test suite, and run `R CMD check` on a clean source package. Document
intentional incompatibilities and migration steps.

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

Use patch releases for compatible fixes. During the 0.x series, use minor
releases for new features and intentional interface changes. Retain older
behavior when the benefit is clear and the implementation remains readable and
inexpensive to maintain or run.

Publish a release only when explicitly requested. A release consists of an
annotated `vMAJOR.MINOR.PATCH` tag and a corresponding hosting-platform release
whose description begins with a concise human-written summary. Verify the tag,
automated checks, and release description before treating publication as
complete.
