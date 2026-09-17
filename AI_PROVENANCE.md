# AI provenance

## Declaration

This repository has received substantial AI coding assistance. AI-generated
output is subject to human selection, review, testing, and approval. Human
contributors remain responsible for the code and documentation they publish.

## Involvement vocabulary

- **Generated:** AI output was accepted substantially as written.
- **Assisted:** AI suggestions were materially revised by a human.
- **Reviewed:** AI evaluated existing work without being its principal
  generator.

## Systems used

| Provider | Tool or agent | Model identifier | Period | Roles | Identification basis |
| --- | --- | --- | --- | --- | --- |
| OpenAI | Codex | GPT-5 | 2026-09 | implementation, testing, documentation, review, refactoring | The active environment identifies the agent as Codex based on GPT-5. |

## Agent registration

An AI agent making a material contribution must read this manifest before
changing the repository. If its provider, tool, and exposed model identifier
are not represented above, it must add a row. If they are represented, it must
extend the period or roles when necessary.

Use only identity information exposed by the environment. Do not infer a model
version or record prompts, private context, conversation logs, or session
identifiers.

Registration describes the system involved. It does not make that system a Git
author or transfer human responsibility. The human contributor should add the
applicable `Assisted-by:` trailer when committing work. Commits and pull
requests made with agent assistance must follow `CONTRIBUTING.md`.

## Human accountability

- A human contributor must understand and approve every submitted change.
- Normal review, testing, security, licensing, and contribution requirements
  apply regardless of how a change was produced.
- AI systems are not Git authors, copyright holders, or DCO signatories.
- `Signed-off-by:` is reserved for a human making the applicable certification.
- An agent must not invent or silently substitute a human, bot, or service
  identity when committing or publishing work.

## Recording limitations

Record the exact provider model identifier when the tool exposes it. Otherwise,
record the visible product or tool label and state that the underlying model is
not disclosed. Do not infer or invent a model version.
