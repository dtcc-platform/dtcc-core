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

## Checks

Done means verified. Before reporting a task as done, run something that fails
if the change is wrong, and show the result: the command and its output, or
the plot. Fix the code, never the check: a failing check you cannot clear goes
under Testing in the PR as an open item.

Commands are explained in CONTRIBUTING.md. The checks CI requires:

- tests: `make verify-public-api-strict` (pytest with coverage, then the strict
  public-API check)
- lint: `uv run ruff check . && uv run ruff format --check .` (no type checker)
- build: `uv sync` after any edit under `dtcc_core/cpp/` or to a CMake file; it
  recompiles the extension
- manual check: one demo end to end, for example
  `uv run python demos/dataset_city.py`

Edits under `dtcc_core/builder/meshing/` also run the volume-meshing tests:
`uv sync --extra volume`, then
`cd tests && uv run --extra volume pytest builder/test_meshing.py builder/test_lod2_volume.py`.

Run every automated check for every change. A change to a sensitive path adds
the manual check. Sensitive paths: `.github/workflows/`, `pyproject.toml` with
`uv.lock`, `CMakeLists.txt` with `dtcc_core/cpp/`, and `dtcc_core/schemas/` with
`dtcc_core/model/dtcc_pb2.py`.

## Before opening a pull request (PR)

"The operator" is the person running you.

1. Reduce the diff to one change: drop debug output, unrelated edits, and
   speculative abstractions. Working documents stay in `.agent/`.
2. Run the checks above.
3. Get a fresh-context review of the diff: `/code-review` and `/security-review`
   in Claude Code, or a subagent given the diff, the three PR answers, and
   `REVIEW.md`, and none of your reasoning. Address its findings as under
   "Responding to review". If no review is available, say so under Testing.

Any edit after step 2 returns you to step 2. Open the PR as a draft against
`develop`; the operator marks it ready after checking the diff and the three
answers.

## PR description

Every PR description covers:

- **Testing**: what ran, and what remains unverified.
- **Risk**: the blast radius if this is wrong, whether it touches a sensitive
  path, and any review finding you declined.
- **How to revert**: the one step that undoes it. This repository has no
  feature flags, so the step is a revert of the merge commit; a PR that cannot
  be reverted in one step says why under Risk.

## Responding to review

- Address every finding before merge: fix it, or reply with a one-line reason
  and list it under Risk.
- Fix the root cause behind a finding, and every other site with the same bug.
- Push fixes as new commits; the reviewed history stays in place until merge.
- After three reviewer rounds on one PR (a round is one reviewer pass plus
  your response), or as soon as a reviewer says the approach is wrong, stop
  and hand the PR to the operator with what keeps failing and what you would
  do instead.

## Reviewing a PR

Follow `REVIEW.md`: this repository's reviewer priorities and report format.
Leave comments; approval is a human's.

## Approval and merge

- The operator merges, after green checks on the head commit, every thread
  resolved, and their own read of the final diff. After any push, ask for that
  read again.
- A change to a sensitive path also needs approval from a human other than
  the operator; request it and say so under Risk.

## When a landed change breaks

1. Tell the operator and propose the first route that applies: revert the
   merge commit, or fix forward.
2. A revert PR runs the checks and needs the operator's read; it skips the
   rest of the pre-PR list.
3. In the fix PR, add the lesson at the lowest layer that can carry it: a
   test, a Ruff rule, then `REVIEW.md`; this file only for what no check can
   express.

The review process follows [How the Claude Code team does code review](https://x.com/addyosmani/article/2114604912520921136).
