# Agent instructions

Use the repository's configured human Git identity for commits. Check
`git var GIT_AUTHOR_IDENT` and `git var GIT_COMMITTER_IDENT` before committing;
both must resolve to the intended human identity from `git config user.name`
and `git config user.email`. Do not override them with an agent identity via
`--author`, `git -c`, Git configuration, or author/committer environment variables.
Never use Claude, Anthropic, Codex, OpenAI, or another agent as author,
committer, or `Co-authored-by`. Do not add generated-by attribution trailers.
If the human identity is missing or ambiguous, ask before committing. Preserve
other contributors' attribution; rewrite existing attribution only when explicitly
requested by the user.

Do not commit implementation plans or temporary working documents, such as
scratch notes, progress logs, review notes, or drafts.

Keep plans in `.agent/plans/` and other working documents in `.agent/work/`.
Both directories are ignored by Git; do not force-add their contents.

Commit lasting project documentation when it is part of the requested change.
Review staged files before committing to catch accidental working documents.

Comments and docstrings describe what the code does now: no authorship notes,
change history, or references to earlier versions, reviews or experiments.

Keep README.md a concise project introduction with documentation links, never
a changelog, work log or design document.
