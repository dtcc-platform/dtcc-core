# Risks

Important findings from merged pull requests that were accepted but left
unfixed. The pull request that fixes one deletes its line. Keep entries sorted
by path.

Entry format; the symbol and PR number are optional:

    - `path/to/file.py` (`symbol`): what is wrong; why it was left open. (#PR)

- `dtcc_core/builder/meshing/tetgen.py` (`build_volume_mesh`): Native city volume builds can drop the input SRS. This predates corner repair and is tracked in [#200](https://github.com/dtcc-platform/dtcc-core/issues/200). (#180)
