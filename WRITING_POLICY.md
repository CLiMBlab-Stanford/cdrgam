# Writing policy

This policy applies when an AI agent creates or revises documentation,
docstrings, or explanatory code comments. Write for researchers and engineers.
Text must be technically accurate, concise, natural, and easy to scan.

## Principles

- Write in plain, direct language.
- State the main point early and explain concepts in the order readers need.
- Prefer short sentences, concrete verbs, and specific nouns.
- Match the project's terminology and documentation conventions.
- Assume the reader is capable but may not know this implementation.
- Include enough context to understand behavior without reading the source.
- Preserve important constraints, qualifications, failure modes, and edge
  cases.
- Use small, valid examples when they clarify behavior faster than prose.
- Keep established mathematical and statistical terms when they are precise.

## Avoid

- Marketing language, hype, cheerleading, and exaggerated claims.
- Jargon that does not improve precision.
- Long noun phrases and unnecessary adjectives or adverbs.
- New labels for concepts that already have standard names.
- Claims that behavior is simple, easy, obvious, or intuitive.
- Stock filler such as “at its core,” “in essence,” “it is important to note,”
  “seamlessly,” “powerful,” or “comprehensive.”
- Contrast formulas such as “not X, but Y” when a direct statement is clearer.
- Excessive headings, lists, parenthetical remarks, em dashes, and rhetorical
  questions.
- Announcing what the documentation will explain.
- Mentioning AI generation in ordinary technical prose. Repository policy and
  commit metadata handle disclosure.

## Documentation requirements

- Describe observable behavior before implementation details.
- Explain purpose, inputs, outputs, side effects, errors, and constraints when
  relevant.
- Use exact names for functions, parameters, classes, commands, and files.
- Do not infer behavior that the code does not establish.
- Focus on intent, contracts, surprising behavior, and decisions callers must
  make instead of restating syntax.
- Keep examples minimal and consistent with the current API.
- Separate distinct ideas into paragraphs.
- Use lists only for genuinely parallel items.

## Docstrings and comments

- Follow R and roxygen2 conventions already used by the package.
- Begin with a brief active-voice summary.
- Document parameters by meaning and constraints, not merely type.
- Describe return values when their meaning is not clear from context.
- Document errors that callers may need to handle.
- Mention mutation, I/O, caching, ordering, units, or lifecycle rules when they
  affect the contract.
- Avoid private implementation details unless they affect callers.
- Keep comments and documentation synchronized with the code.

## Review method

1. Read the relevant code and nearby documentation.
2. Identify the audience and required facts.
3. Draft the smallest complete explanation.
4. Check every technical claim against code or supplied source material.
5. Edit for clear structure, brevity, and natural rhythm.
6. Remove generic filler.
7. Confirm that terminology and formatting match the repository.
