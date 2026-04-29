# Merge `develop` into `feature/lod2-reconstruction` — Trial-Branch Runbook

**Date:** 2026-04-29
**Goal:** Bring 98 commits of `develop` into the LoD2 feature branch, but evaluate the cost on a throwaway branch first. If it works, fast-forward the real branch. If it doesn't, throw the trial away — the real branch stays untouched.

**Pre-flight expectation:** 3 files conflict (`.gitignore`, `pyproject.toml`, `dtcc_core/model/mixins/city/builder_mixin.py`). All three are tractable.

---

## 0. Pre-flight checks

Run from the repo root: `/Users/vasnas/scratch/lod3/temp/dtcc-core`.

```bash
# Must be on the feature branch with a clean tree
git status
```

Expected: `On branch feature/lod2-reconstruction` and `nothing to commit, working tree clean`.

If untracked files exist (e.g. the brainstorm `.txt`), that's fine — they don't affect the merge.

```bash
# Tag the current state as a safety anchor before doing anything
git tag pre-develop-merge-2026-04-29
```

Now if anything goes catastrophically wrong, `git reset --hard pre-develop-merge-2026-04-29` puts you back where you are right now.

---

## 1. Create the trial branch

```bash
git checkout -b feature/lod2-reconstruction-merge-trial
```

Confirm:

```bash
git status
git log --oneline -1
```

Expected: branch is `feature/lod2-reconstruction-merge-trial`, HEAD same commit as the original feature branch.

---

## 2. Fetch develop and attempt the merge

```bash
git fetch origin develop
git merge origin/develop
```

**Three outcomes possible. Find yours below.**

---

## 3a. Outcome A — Clean merge (git auto-resolves everything)

You'll see a successful merge message ending with `Merge made by the 'ort' strategy.` or similar. No "CONFLICT" lines.

Skip to **Section 4 — Verify**.

---

## 3b. Outcome B — Conflicts (most likely; expected on the 3 files)

You'll see lines like:

```
Auto-merging .gitignore
CONFLICT (content): Merge conflict in .gitignore
Auto-merging pyproject.toml
CONFLICT (content): Merge conflict in pyproject.toml
Auto-merging dtcc_core/model/mixins/city/builder_mixin.py
CONFLICT (content): Merge conflict in dtcc_core/model/mixins/city/builder_mixin.py
Automatic merge failed; fix conflicts and then commit the result.
```

Resolve each. Guidance per file follows.

### 3b.1 — `.gitignore`

Open it. You'll see conflict markers around the bottom of the file:

```
<<<<<<< HEAD
sandbox/snapping

# Pipeline run artefacts (OBJ/VTK/VTU/JSON/HTML/etc produced by sandbox scripts)
sandbox/output/

# uv lockfile (library convention — reproduce via pyproject, not a pinned tree)
uv.lock
=======
sandbox/snapping

# Generated experiment artifacts
sandbox/output/
sandbox/output_3d_*/
sandbox/output_stage_dive_artifacts_*/
sandbox/tetgen_fail*
sandbox/tetgen-tmpfile_skipped*
tetgen_fail*
tetgen-tmpfile_skipped*
benchmarks/__pycache__/
benchmarks/output*/
>>>>>>> origin/develop
```

Resolution: keep both sets, dedupe `sandbox/output/`. Replace the entire conflict block with:

```
sandbox/snapping

# Generated experiment artifacts
sandbox/output/
sandbox/output_3d_*/
sandbox/output_stage_dive_artifacts_*/
sandbox/tetgen_fail*
sandbox/tetgen-tmpfile_skipped*
tetgen_fail*
tetgen-tmpfile_skipped*
benchmarks/__pycache__/
benchmarks/output*/

# uv lockfile (library convention — reproduce via pyproject, not a pinned tree)
uv.lock
```

Develop also added `temp/` higher up in the file — accept their version of that line (it's not a conflict zone, but check there's no extra HEAD/======/>>>>>> markers anywhere else in the file).

Mark resolved:

```bash
git add .gitignore
```

### 3b.2 — `pyproject.toml`

Open it. Conflicts will be in the `dependencies = [...]` block and `[project.optional-dependencies]` block.

**Strategy:** keep develop's restructure (it's bigger and more consequential), then re-apply our `open3d` line.

Specifically:
- In `dependencies`, develop removed `dtcc-pyspade-native` and added `dtcc-mesher`. Keep their version. Then add `"open3d>=0.17.0",` somewhere in the list (alphabetical-ish placement matches their style).
- In `[project.optional-dependencies]`, develop has:
  ```toml
  test = ["pytest", "pytest-cov", "httpx>=0.27"]
  spade = ["dtcc-pyspade-native == 0.1.3"]
  remote = ["httpx>=0.27"]
  ```
  Keep their version verbatim. We don't need to add an `[ml]` extras here (Track 2 spec calls for it but Track 2 hasn't started; add later when the classifier lands).
- Develop bumped `version = "0.9.6dev"` → `"0.9.8dev"`. Keep theirs.
- Develop changed `[build-system].requires` to drop `dtcc-pyspade-native` and add `ninja`. Keep theirs.

Verify the file has no remaining `<<<<<<<`, `=======`, or `>>>>>>>` markers:

```bash
grep -n '<<<<<<<\|=======\|>>>>>>>' pyproject.toml
```

Expected: empty output.

Mark resolved:

```bash
git add pyproject.toml
```

### 3b.3 — `dtcc_core/model/mixins/city/builder_mixin.py`

This one should auto-merge or be a clean structural conflict. Our change inserted a brand-new method `build_lod2_buildings`; develop modified existing methods. The conflict, if any, is just about where in the file the new method lands.

Open the file and look for conflict markers:

```bash
grep -n '<<<<<<<\|=======\|>>>>>>>' dtcc_core/model/mixins/city/builder_mixin.py
```

If empty: git auto-merged. Move on.

If markers exist: ensure the `build_lod2_buildings` method we added is preserved (search for `def build_lod2_buildings`), and develop's parameter additions to `build_surface_mesh`, `build_flat_mesh`, `build_volume_mesh`, `build_terrain` (look for `mesher`, `pipeline_mode`, `show_footprints`, `top_cap_max_mesh_size`, `max_volume`) are also preserved. Both should coexist — they touch different methods.

Verify clean:

```bash
grep -n '<<<<<<<\|=======\|>>>>>>>' dtcc_core/model/mixins/city/builder_mixin.py
```

Mark resolved:

```bash
git add dtcc_core/model/mixins/city/builder_mixin.py
```

### 3b.4 — Finalize the merge commit

```bash
git status
```

Expected: `All conflicts fixed but you are still merging.`

```bash
git merge --continue
```

This drops you into your editor with a default merge commit message. Accept it (`:wq` in vi, `Ctrl+X` then `Y` in nano) or amend.

Confirm:

```bash
git log --oneline -3
```

Expected top commit: `Merge remote-tracking branch 'origin/develop' into feature/lod2-reconstruction-merge-trial`.

---

## 3c. Outcome C — Merge looks much worse than expected

If you see conflicts in many more than 3 files, or in files we never touched, abort:

```bash
git merge --abort
```

This puts you back to where you were before `git merge origin/develop` started.

Then go to **Section 6 — Abandon trial**.

---

## 4. Verify the merge result

### 4.1 — Refresh dependencies

Develop dropped `dtcc-pyspade-native` from main deps and added `dtcc-mesher`. The venv must follow.

```bash
.venv/bin/pip install -e . 2>&1 | tail -10
```

Watch for errors. If you actually use `pyspade` in any sandbox script you call, also install the spade extras:

```bash
.venv/bin/pip install -e .[spade]
```

### 4.2 — Run the LoD2-specific tests

```bash
.venv/bin/pytest \
    tests/builder/test_lod2_pipeline.py \
    tests/builder/test_lod1_regression.py \
    tests/builder/test_roof_config.py \
    tests/builder/test_roof_detection.py \
    tests/builder/test_roof_detection_timing.py \
    tests/builder/test_roof_classification.py \
    tests/builder/test_roof_geometry.py \
    tests/builder/test_shell_validation.py \
    tests/builder/test_extract_roof_points_fix.py \
    tests/builder/test_build_lod2_diagnostics.py \
    tests/builder/evaluation/ \
    tests/model/test_multisurface_semantics.py \
    tests/model/test_pointcloud_normals.py \
    -v
```

Expected: all pass. If any fail, look at the failure — most likely cause is develop changing something we depend on (e.g., a pointcloud or surface API).

### 4.3 — Run the full suite

```bash
.venv/bin/pytest tests/ 2>&1 | tail -3
```

Expected: counts higher than 689 (develop added ~30 test files), zero failures.

If develop's tests fail (not ours), that's their problem to fix — not a merge issue. But it's still worth investigating.

### 4.4 — Smoke-run the eval harness CLI

```bash
.venv/bin/python sandbox/evaluate_lod2.py \
    --dataset tests/builder/evaluation/fixtures/minimal_dataset \
    --out /tmp/post_merge_check
```

Expected: `Evaluated 1 buildings` and CSV/JSON written under `/tmp/post_merge_check/`.

---

## 5. Decision: keep or discard

### 5a. If verification passed — fast-forward the real branch

```bash
git checkout feature/lod2-reconstruction
git merge --ff-only feature/lod2-reconstruction-merge-trial
```

If `--ff-only` succeeds, the real branch now has the merged history. Confirm:

```bash
git log --oneline -3
```

Expected top commit: the merge commit from the trial.

Clean up the trial branch:

```bash
git branch -d feature/lod2-reconstruction-merge-trial
```

Push if you want it on origin:

```bash
git push origin feature/lod2-reconstruction
```

(No force-push needed — fast-forward only.)

### 5b. If verification failed — discard

```bash
git checkout feature/lod2-reconstruction
git branch -D feature/lod2-reconstruction-merge-trial
```

Real branch is exactly where it was. The `pre-develop-merge-2026-04-29` tag is no longer needed but harmless; delete if you want:

```bash
git tag -d pre-develop-merge-2026-04-29
```

Decide later whether to retry the merge after develop stabilizes, do a piecewise pick of just the parts you need, or stay on the current divergence.

---

## 6. Abandon trial (used by Outcome C above)

Already on the trial branch with the merge aborted (or never started):

```bash
git checkout feature/lod2-reconstruction
git branch -D feature/lod2-reconstruction-merge-trial
```

Tag stays for safety; remove with `git tag -d pre-develop-merge-2026-04-29` if desired.

---

## Reference: which files to expect conflicts in

```
.gitignore
pyproject.toml
dtcc_core/model/mixins/city/builder_mixin.py
```

If conflict-marker search reveals others:

```bash
grep -rn '<<<<<<<\|=======\|>>>>>>>' --include='*.py' --include='*.toml' --include='*.md' --exclude-dir=.git --exclude-dir=.venv 2>&1 | head -20
```

Investigate before continuing — that's a signal something unexpected is going on.

---

## Reference: emergency rollback

If the trial branch got fast-forwarded into the real branch and tests immediately reveal a problem you can't fix:

```bash
git checkout feature/lod2-reconstruction
git reset --hard pre-develop-merge-2026-04-29
git push --force-with-lease origin feature/lod2-reconstruction
```

This is destructive — only use if you're sure. The `--force-with-lease` flag is safer than `--force` because it refuses if origin moved between your fetch and push.
