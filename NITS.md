# Nits

Nit findings from merged pull requests that were accepted but left unfixed.
The pull request that fixes one deletes its line. Keep entries sorted by path.

Entry format; the symbol and PR number are optional:

    - `path/to/file.py` (`symbol`): what is wrong; why it was left open. (#PR)

- `dtcc_core/datasets/transport/vasttrafik.py` (`_fetch_token`, `_fetch_positions`): a failing live run under pytest's default or `-l` traceback prints these frames' credential arguments; the live workflow runs with `--tb=short`, which does not, so it was left open. (#139)
