# Simplify dataset plotting modes and add rich smoke preview

Status: completed
Created: 2026-07-01
Suggested path: `.agent/plans/2026-07-01-dataset-plot-modes-and-smoke-preview.md`

This plan is self-contained enough for Codex to implement without needing the original discussion. Keep the implementation practical, explicit, and scoped. The task is to simplify the public plotting API for DTCC datasets and make the smoke dataset the first polished reference implementation.

## Goal

Implement a simple three-mode plotting contract for dataset plotting in `dtcc-core`, starting with `datasets.smoke.plot(...)`:

1. `mode="artifact"`: bare renderable artifact for tangible-table export/publishing. No baked-in narrative, titles, axes, or frontend-owned UX.
2. `mode="preview"`: rich Python preview of the intended tangible-table/museum-style experience, including narrative, legend, annotations, key facts, and limitations.
3. `mode="plot"`: ordinary Python/Matplotlib plot with title, axes, and simple legend/colorbar.

The default developer experience should become:

```python
dtcc.datasets.smoke.plot(bounds=bounds)
```

which should produce a polished preview figure, not a plain technical plot.

At the same time, preserve the existing practical embedding behavior:

```python
fig, axes = plt.subplots(1, 2)
dtcc.datasets.smoke.plot(bounds=bounds, ax=axes[0], product="slice", show=False)
```

When `ax` is supplied and `mode` is omitted, default to `mode="plot"` so existing Matplotlib subplot usage remains intuitive.

The smoke dataset should become the reference “all bells and whistles” synthetic dataset for Dataset v2 presentation: deterministic, visually rich, metadata-rich, and useful as a local preview of what the tangible table may eventually show.

## Non-goals

- Do not implement the tangible-table frontend.
- Do not hardcode final table UX decisions that should belong to table designers/developers.
- Do not bake narrative panels, legends, or metadata into artifacts intended for table consumption.
- Do not redesign all dataset plotting in `dtcc-core`; implement the smoke reference path first and keep reusable helpers modest.
- Do not complete Dataset v2 packaging/export/publish phases beyond what is necessary for plotting/export preview consistency.
- Do not introduce new heavy dependencies, web assets, custom fonts, GUI frameworks, or non-Matplotlib rendering stacks.
- Do not remove the existing `profile` argument immediately; keep it as a backward-compatible alias during this change.
- Do not change the mathematical smoke field unless needed for visual presentation defaults.

## Background

The current smoke dataset already has a useful technical base. It exposes products such as `field`, `slice`, and `streamlines`; supports formats including `pb`, `vtu`, `geojson`, `png`, and `mp4`; and has render options for table/Python profiles, theme, colormap, dimensions, legends, titles, line glow, and animation.

The current `SmokeDataset.plot()` is a thin Matplotlib wrapper. It defaults to a slice product, sets `profile="python"`, enables a legend, assigns a simple title, builds a visualization product, and delegates to the generic product renderer. This is useful for technical inspection, but it feels like a normal scientific plot rather than a preview of an audience-facing tangible-table experience.

The current `profile` split is also not the right public abstraction. The user-facing decision is not really “profile plus view”; it is the intended output mode:

- artifact for table consumption;
- preview for Python-side experience preview;
- plot for ordinary Python plotting.

Dataset v2 already distinguishes metadata, provenance, and presentation. Presentation is the audience/UX layer and is intended to support table cards, narrative panels, Python plotting, public demos, legends, annotations, and view hints. The smoke preview should consume this kind of presentation information rather than hardcoding unexplained text in plotting logic.

Important decisions already made:

- Use one public option: `mode`.
- Supported modes are exactly `"artifact"`, `"preview"`, and `"plot"`.
- Keep `product` as the existing data/product selection concept, but give preview/artifact a good default visual composition for smoke.
- Keep `profile` temporarily as an alias:
  - `profile="table"` maps to `mode="artifact"` when `mode` is omitted.
  - `profile="python"` maps to `mode="plot"` when `mode` is omitted.
- If both `mode` and `profile` are supplied and conflict, fail loudly.
- Default `datasets.smoke.plot(bounds=bounds)` should be `mode="preview"`.
- Default `datasets.smoke.plot(..., ax=ax)` should be `mode="plot"` to preserve subplot ergonomics.

## Acceptance criteria

- [x] `dtcc.datasets.smoke.plot(bounds=bounds, show=False)` produces the rich preview mode by default and returns a Matplotlib `Figure`.
- [x] `dtcc.datasets.smoke.plot(bounds=bounds, mode="preview", show=False)` produces the same rich preview mode and returns a Matplotlib `Figure`.
- [x] `dtcc.datasets.smoke.plot(bounds=bounds, mode="plot", show=False)` produces an ordinary Python plot and returns a Matplotlib `Axes`.
- [x] `dtcc.datasets.smoke.plot(bounds=bounds, mode="artifact", show=False)` produces a bare artifact-style plot and returns a Matplotlib `Axes`.
- [x] `dtcc.datasets.smoke.plot(bounds=bounds, ax=ax, show=False)` defaults to `mode="plot"` and returns the supplied `Axes`.
- [x] `dtcc.datasets.smoke.plot(bounds=bounds, mode="preview", ax=ax, show=False)` fails with a clear `ValueError` explaining that preview mode owns the full figure layout and cannot render into a single axes.
- [x] `profile="table"` remains accepted for smoke plotting when `mode` is omitted and maps to `mode="artifact"`.
- [x] `profile="python"` remains accepted for smoke plotting when `mode` is omitted and maps to `mode="plot"`.
- [x] Supplying conflicting `mode` and `profile` values fails with a clear `ValueError`; no silent precedence rule is allowed.
- [x] Invalid `mode` values fail with a clear error listing `artifact`, `preview`, and `plot`.
- [x] Preview mode includes: main visual, headline, short summary, legend, at least two narrative cards, at least two meaningful annotations/callouts, key facts, and limitations/footer.
- [x] Artifact mode includes no title, no axes, no Matplotlib colorbar, no narrative panel, no callout text, and no frontend-owned layout.
- [x] Plot mode preserves ordinary Python plot behavior: title, axis labels, optional simple legend/colorbar, and compatibility with supplied `ax`.
- [x] Smoke dataset `describe()` exposes structured presentation metadata for the preview.
- [x] Dataset-produced smoke objects carry presentation metadata through `dataset_context` where the generic Dataset v2 context mechanism supports it.
- [x] Existing smoke PNG/MP4 artifact export behavior remains compatible unless explicitly changed by the plan.
- [x] Tests cover all three modes, invalid/conflicting mode/profile inputs, and `ax` behavior.
- [x] `demos/smoke.py` is updated to demonstrate the three modes without requiring manual code edits.
- [x] No new external runtime dependency is added.

## Fail-loud requirements

- Required item: `bounds`
  - Valid when: a valid 4-value or 6-value DTCC bounds sequence/object satisfying existing `DatasetBaseArgs` validation.
  - Invalid/missing behavior: fail with the existing clear bounds validation error.
  - Silent fallback forbidden: yes.

- Required item: `mode`
  - Valid when: omitted, `"artifact"`, `"preview"`, or `"plot"`.
  - Invalid/missing behavior: omitted is resolved by the rules below; invalid explicit values fail with a clear `ValueError` listing valid modes.
  - Silent fallback forbidden: yes for invalid explicit values.

- Required item: mode resolution when `mode` is omitted
  - Valid when:
    - `ax is None` and no legacy `profile` was supplied: resolve to `"preview"`.
    - `ax is not None` and no legacy `profile` was supplied: resolve to `"plot"`.
    - `profile="table"`: resolve to `"artifact"`.
    - `profile="python"`: resolve to `"plot"`.
  - Invalid/missing behavior: any unsupported profile value fails via existing validation; unsupported combinations fail clearly.
  - Silent fallback forbidden: yes.

- Required item: `mode` and `profile` compatibility
  - Valid when:
    - `mode="artifact"` with `profile="table"`.
    - `mode="plot"` with `profile="python"`.
    - `mode="preview"` with no `profile`.
  - Invalid/missing behavior: conflicting combinations fail with `ValueError`, for example `mode="preview", profile="table"` or `mode="artifact", profile="python"`.
  - Silent fallback forbidden: yes.

- Required item: preview mode layout ownership
  - Valid when: `mode="preview"` and `ax is None`.
  - Invalid/missing behavior: `mode="preview"` with `ax` fails with `ValueError` explaining that preview mode creates a full figure layout.
  - Silent fallback forbidden: yes.

- Required item: Matplotlib
  - Valid when: Matplotlib can be imported by the existing lazy import mechanism.
  - Invalid/missing behavior: fail with the existing clear Matplotlib-required error.
  - Silent fallback forbidden: yes.

- Required item: smoke presentation metadata
  - Valid when: smoke descriptor provides non-empty headline, summary, narrative, legend, annotations, key points, and limitations for preview mode.
  - Invalid/missing behavior: smoke preview should fail during development/tests if required smoke presentation fields are missing or malformed. For future datasets, the generic preview renderer may degrade only when explicitly designed to do so.
  - Silent fallback forbidden: yes for smoke reference data.

- Required item: annotation coordinates in smoke presentation
  - Valid when: relative annotation positions are two floats in `[0.0, 1.0]`.
  - Invalid/missing behavior: fail with `ValueError` in preview mode rather than drawing misplaced callouts.
  - Silent fallback forbidden: yes for smoke reference annotations.

- Required item: artifact mode purity
  - Valid when: artifact mode only renders data layers and frontend-safe visual styling.
  - Invalid/missing behavior: tests should fail if artifact mode includes title, axes, narrative text, callout text, or metadata panels.
  - Silent fallback forbidden: yes.

## CLI ergonomics requirements

Not applicable.

This task changes Python plotting behavior and does not create or modify a human-facing CLI tool.

## Relevant files

Likely files to inspect or modify:

- `dtcc_core/datasets/smoke.py`: primary implementation target for smoke plot modes, mode resolution, smoke presentation metadata, and smoke-specific preview composition.
- `dtcc_core/plotting/renderers.py`: current generic raster/product rendering. Reuse existing `draw_product`, `_draw_slice`, and `_draw_streamlines` where practical. Avoid breaking existing render paths.
- `dtcc_core/plotting/options.py`: current raster/video render options. Add only minimal mode-related fields if they materially simplify implementation; otherwise keep mode handling in smoke plotting helpers.
- `dtcc_core/plotting/products.py`: current `SliceProduct` and `StreamlineProduct`. Add a reusable composite/overlay product only if it is cleaner than smoke-local composition.
- `dtcc_core/plotting/style.py`: current DTCC colors, themes, colorbar styling, metadata boxes, and helpers. Reuse palette and theme values.
- `dtcc_core/datasets/dataset.py`: generic dataset context creation. Consider consuming descriptor-level `presentation` data here so smoke-produced objects expose richer `obj.presentation`.
- `dtcc_core/datasets/schema.py`: DatasetPresentation schema already exists. Modify only if needed; prefer using existing fields.
- `demos/smoke.py`: update the demo to show artifact, preview, and plot modes.
- `DESIGN.md`: reference only; do not rewrite unless a small docs update is needed.
- `tests/`: discover existing dataset or plotting tests. Add or update focused tests for smoke plot modes. If no suitable file exists, add `tests/datasets/test_smoke_plot_modes.py` or another path consistent with repository conventions.

If relevant tests are not obvious, discover them with:

```bash
find tests -iname '*smoke*' -o -iname '*plot*' -o -iname '*dataset*'
```

## Implementation approach

### Public API design

Add a single public plotting option to smoke plotting:

```python
SmokePlotMode = Literal["artifact", "preview", "plot"]

def plot(self, ax=None, show: bool = True, mode: SmokePlotMode | None = None, **kwargs):
    ...
```

Do not expose both `profile` and `view` as public concepts for plotting. `profile` may remain in `SmokeArgs` and render options for backward compatibility and export internals, but smoke plotting should treat it as a legacy alias.

Mode semantics:

```text
artifact = bare table/export artifact preview
preview  = rich Python approximation of tangible-table/museum UX
plot     = ordinary Python/Matplotlib plot
```

Mode defaulting:

```python
if mode is None and legacy_profile is None:
    mode = "plot" if ax is not None else "preview"
```

Legacy profile mapping:

```python
profile="table"  -> mode="artifact"
profile="python" -> mode="plot"
```

Conflict handling:

```python
mode="artifact" is compatible only with profile omitted or profile="table"
mode="plot"     is compatible only with profile omitted or profile="python"
mode="preview"  is compatible only with profile omitted
```

Prefer a small helper:

```python
def _resolve_plot_mode(*, mode, profile, ax) -> SmokePlotMode:
    ...
```

and a second helper:

```python
def _render_options_for_mode(args: SmokeArgs, mode: SmokePlotMode) -> RasterRenderOptions:
    ...
```

### Product composition

For smoke, the strongest artifact/preview visual should be an overview composition:

```text
speed slice background + velocity/speed streamlines overlay
```

Do not force the user to understand this composition in the common case. If `product` is omitted:

- `mode="preview"`: use the overview composition.
- `mode="artifact"`: use the overview composition.
- `mode="plot"`: preserve current behavior and default to `product="slice"`.

If `product="slice"` or `product="streamlines"` is explicitly supplied:

- `mode="plot"`: render the requested single product using current behavior.
- `mode="artifact"`: render the requested product without axes/title/legend/narrative.
- `mode="preview"`: either:
  - use the requested product as the main visual if this is straightforward; or
  - fail with a clear error if preview currently supports only the overview composition.

Preferred first implementation: support `product="slice"`, `product="streamlines"`, and omitted/overview in preview mode, but keep the omitted overview as the polished path.

Avoid adding `product="overview"` to the dataset build API unless Codex determines it can be done cleanly without disrupting existing serialization semantics. The overview is primarily a visualization composition, not necessarily a new serialized dataset product.

### Preview styling requirements

Preview mode should be visually much better than the current plot. Treat it as a compact museum/kiosk exhibit panel, implemented with Matplotlib primitives and existing DTCC colors.

Use the existing DTCC palette as the visual foundation:

```text
background: #0F0F14
main visual background: #101016
panel: #1B1B22
card/panel secondary: #24242B
text: #FFFFFF
muted text: #C9C9D1
grid/edge: #3A3A42
yellow accent: #FADA36
teal accent: #78C8BE
orange accent: #E35A1D
```

Typography:

- Use Matplotlib/default installed fonts only, preferably default sans-serif/DejaVu Sans.
- Do not add custom font files.
- Use strong typographic hierarchy:
  - headline: large and bold;
  - section/card headings: medium bold;
  - body text: smaller, muted, readable;
  - footer: small and quiet.
- Avoid dense text. Use short sentences and wrap text manually where necessary.

Recommended preview layout for 16:9 figure:

```text
Canvas: dark full figure, no global axes.
Left visual card: roughly 68% width, 84% height, rounded panel, subtle shadow.
Right narrative panel: roughly 23-28% width, 84% height, rounded panel.
Bottom footer: quiet provenance/limitation line.
```

Concrete relative layout guidance:

```text
main visual card: x=0.04, y=0.08, w=0.66, h=0.84
right panel:       x=0.735, y=0.08, w=0.225, h=0.84
footer:            y≈0.03
```

These exact values may be adapted to cleaner Matplotlib code, but the final composition should have the same visual balance: immersive data on the left, concise story on the right.

Main visual styling:

- Draw the scalar speed slice as a full-bleed image using the DTCC numeric colormap.
- Overlay streamlines using a two-pass glow:
  - wide low-alpha yellow/teal glow underneath;
  - narrow high-alpha bright line on top.
- Use rounded-card framing and subtle shadow around the visual area.
- Add a very subtle grid/fiducial overlay only if it improves the table-preview feel; it must not look like ordinary plot axes.
- Do not show standard axis ticks in preview mode.
- Prefer an integrated horizontal color ramp instead of a standard Matplotlib colorbar.
- Keep the color ramp compact, inside or just below the main visual card, with labels such as `Smoke speed`, `slow`, `fast`, and `[m/s]`.

Annotation/callout styling:

- Use at most 2-3 callouts in the default smoke preview.
- Use small accent dots/circles anchored to relative coordinates over the main visual.
- Use thin leader lines or arrows from the card label to the anchor.
- Use short labels, for example:
  - `Fast corridor`
  - `Recirculation`
  - `Slice plane`
- Use short descriptions, for example:
  - `warmer color + tighter flow`
  - `looping trace suggests mixing`
- Callouts should be visible but not dominate the data.
- Validate callout coordinates and fail loudly for malformed smoke defaults.

Right panel content:

1. Headline:

   ```text
   Synthetic Urban Smoke Flow
   ```

2. Summary:

   ```text
   A deterministic smoke-test dataset showing how a velocity field can be published, previewed, animated, and consumed by the tangible table.
   ```

3. Chips/badges:

   ```text
   simulation
   EPSG:3006
   synthetic
   loopable
   ```

4. Narrative cards:

   ```text
   What you are seeing
   The background shows smoke speed on a slice through the volume. The bright curves trace the local flow direction.

   How to read it
   Warm colors indicate faster motion. Curved and closed traces suggest recirculation and mixing zones.

   Why it matters
   This is a safe synthetic stand-in for environmental simulation data while the table interaction model is developed.
   ```

5. Key facts/stat strip:

   Examples:

   ```text
   64 flow lines
   8s loop
   PNG preview
   MP4 motion
   ```

   Compute facts from actual args/products where possible. Do not fake numbers. If a fact cannot be computed, omit that fact rather than using placeholder text.

6. Footer/limitations:

   ```text
   synthetic demonstration data · not a validated CFD simulation · not a forecast
   ```

Preview background and panels:

- Use rounded rectangles for visual and text panels.
- Use subtle edge colors from the DTCC grid color.
- Use slightly transparent dark panels only if this renders reliably in Matplotlib backends.
- Avoid clutter. The preview should feel like an exhibit card, not a dashboard.

Preview figure size:

- Default to 16:9, e.g. `figsize=(16, 9)`.
- Respect existing `width`, `height`, and `dpi` options if they are passed and can be cleanly mapped to `figsize`.
- If the user passes obviously invalid dimensions, rely on existing validation.

### Artifact mode styling

Artifact mode is the visual layer intended for the tangible-table package. It must be frontend-safe.

Default artifact should be:

```text
speed slice + streamline overlay
```

unless the user explicitly requests `product="slice"` or `product="streamlines"`.

Artifact mode requirements:

- No title.
- No axes.
- No tick labels.
- No standard colorbar.
- No narrative panel.
- No annotations/callouts with text.
- No metadata footer.
- Full-bleed or projection-safe fill.
- Use dark background by default.
- Use DTCC numeric colormap for scalar layer.
- Use glow streamlines where visually appropriate.
- Preserve existing `width`, `height`, `dpi`, `theme`, `background`, `transparent`, `cmap`, `vmin`, `vmax`, `line_width`, `line_alpha`, and `glow` behavior where possible.

### Plot mode styling

Plot mode should remain ordinary, predictable Matplotlib behavior.

Requirements:

- Use current product renderer where possible.
- Title defaults to `DTCC Smoke Slice` or `DTCC Smoke Streamlines` depending on product.
- Axis labels are visible.
- Simple Matplotlib colorbar/legend is allowed.
- Respect supplied `ax`.
- Return `Axes`.
- Preserve current default product `slice` unless explicitly changed.

### Presentation metadata

Add smoke presentation metadata in one place, preferably a helper in `smoke.py`:

```python
def _smoke_presentation() -> dict[str, Any]:
    return {
        "headline": "Synthetic Urban Smoke Flow",
        "summary": ...,
        "narrative": [...],
        "key_points": [...],
        "legend": {...},
        "annotations": [...],
        "view_hints": {...},
        "limitations": [...],
    }
```

`SmokeDataset.describe()` should include this under `metadata["presentation"]` or top-level `"presentation"`, whichever best fits existing Dataset v2 conventions in the repository. Prefer top-level `"presentation"` because `DatasetPresentation` is distinct from factual metadata.

If feasible, update `DatasetDescriptor.create_context()` so descriptor-level presentation fields populate `DatasetContext.presentation` instead of always using only the title and description. Keep this generic and conservative:

- If descriptor has no presentation payload, current behavior remains unchanged.
- If descriptor has presentation payload, merge it with default headline/summary.
- Invalid presentation payload should fail clearly during context creation for the dataset that provided it.

Suggested smoke presentation payload:

```python
{
    "headline": "Synthetic Urban Smoke Flow",
    "summary": (
        "A deterministic smoke-test dataset showing how a velocity field can be "
        "published, previewed, animated, and consumed by the tangible table."
    ),
    "narrative": [
        {
            "heading": "What you are seeing",
            "body": (
                "The background shows smoke speed on a slice through the volume. "
                "The bright curves trace the local flow direction."
            ),
        },
        {
            "heading": "How to read it",
            "body": (
                "Warm colors indicate faster motion. Curved and closed traces "
                "suggest recirculation and mixing zones."
            ),
        },
        {
            "heading": "Why it matters",
            "body": (
                "This is a safe synthetic stand-in for environmental simulation "
                "data while the table interaction model is developed."
            ),
        },
    ],
    "key_points": [
        "Deterministic synthetic vector field",
        "Supports field, slice, and streamline products",
        "Exports GeoJSON, PNG, and MP4 artifacts",
        "Designed to exercise Dataset v2 presentation metadata",
    ],
    "legend": {
        "title": "Smoke speed",
        "unit": "m/s",
        "colormap": "dtcc",
        "overlays": [
            {"label": "Streamlines", "meaning": "local flow direction"},
            {"label": "Glow", "meaning": "table-friendly motion emphasis"},
        ],
    },
    "annotations": [
        {
            "label": "Fast corridor",
            "description": "Warmer color and tighter flow indicate faster motion.",
            "position": [0.58, 0.64],
            "coordinates": "relative",
        },
        {
            "label": "Recirculation",
            "description": "Curving traces suggest local mixing.",
            "position": [0.31, 0.38],
            "coordinates": "relative",
        },
    ],
    "view_hints": {
        "default_plot_mode": "preview",
        "artifact_mode": {
            "default_visual": "speed_slice_with_streamlines",
            "include_text": False,
            "include_axes": False,
        },
        "preview_mode": {
            "layout": "exhibit_panel",
            "aspect": "16:9",
            "main_visual": "speed_slice_with_streamlines",
        },
        "plot_mode": {
            "default_product": "slice",
            "include_axes": True,
            "include_simple_legend": True,
        },
    },
    "limitations": [
        "Synthetic demonstration data; not a validated CFD simulation.",
        "Velocity values are illustrative and should not be interpreted as a forecast.",
    ],
}
```

### Validation and error handling strategy

- Validate `mode` before calling the existing dataset argument model.
- Strip or resolve legacy `profile` before building render options.
- Do not allow `preview` to render into an existing `ax`.
- Do not silently ignore conflicting `mode`/`profile` combinations.
- Do not silently produce a blank preview if presentation metadata is malformed.
- For future generic preview support, graceful degradation may be acceptable, but smoke is the reference dataset and should fail loudly during development.

### Compatibility concerns

- Keep `profile` in `SmokeArgs` for now.
- Do not break existing export calls that pass `profile="table"` or `profile="python"`.
- Existing code that calls `datasets.smoke.plot(..., ax=ax)` should continue to produce a normal plot unless it explicitly asks for `mode="preview"`.
- Existing code that calls `datasets.smoke.plot(bounds=bounds)` without `ax` will intentionally get a richer default preview. This is the desired behavior change.
- Existing `render_product_png` and `render_product_mp4` should remain usable for `slice` and `streamlines`.

### Testing strategy

Use non-interactive Matplotlib backend in tests if the repository does not already configure one:

```python
import matplotlib
matplotlib.use("Agg")
```

Test behavior, not pixel-perfect images.

Suggested tests:

- `test_smoke_plot_default_preview_returns_figure`
- `test_smoke_plot_preview_rejects_ax`
- `test_smoke_plot_with_ax_defaults_to_plot_mode`
- `test_smoke_plot_mode_plot_returns_axes`
- `test_smoke_plot_mode_artifact_returns_axes_without_axis`
- `test_smoke_plot_invalid_mode_fails`
- `test_smoke_plot_profile_alias_table_maps_to_artifact`
- `test_smoke_plot_profile_alias_python_maps_to_plot`
- `test_smoke_plot_conflicting_mode_profile_fails`
- `test_smoke_describe_includes_presentation`
- `test_smoke_dataset_context_includes_presentation_if_context_is_supported`

Avoid image snapshot tests unless the repository already has a robust image comparison setup.

## Milestones

### Milestone 1: Add mode resolver and smoke presentation metadata

Expected changes:

- Add `SmokePlotMode = Literal["artifact", "preview", "plot"]` in `dtcc_core/datasets/smoke.py`.
- Add `_resolve_plot_mode(...)` helper.
- Add `_smoke_presentation()` helper.
- Extend `SmokeDataset.describe()` to include structured presentation metadata.
- If feasible, update `DatasetDescriptor.create_context()` to consume descriptor-level presentation data.
- Add focused tests for mode resolution, invalid/conflicting `mode`/`profile`, and smoke presentation metadata.

Verification:

```bash
python -m pytest tests -k "smoke and presentation"
python -m pytest tests -k "smoke and mode"
```

Status: completed

### Milestone 2: Implement artifact and plot modes

Expected changes:

- Refactor current `SmokeDataset.plot()` into explicit mode branches.
- Implement `mode="plot"` using the current product renderer behavior.
- Implement `mode="artifact"` as a bare frontend-safe visual:
  - no axes;
  - no title;
  - no legend/colorbar;
  - no narrative;
  - no annotation text.
- Implement a smoke overview visual composition for artifact mode when `product` is omitted:
  - speed slice background;
  - streamline overlay;
  - glow enabled by default.
- Preserve support for explicit `product="slice"` and `product="streamlines"`.
- Keep `ax` support for `artifact` and `plot` modes.

Verification:

```bash
python -m pytest tests -k "smoke and plot"
```

Manual check:

```bash
python - <<'PY'
import matplotlib
matplotlib.use("Agg")
import dtcc_core as dtcc
bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
ax = dtcc.datasets.smoke.plot(bounds=bounds, mode="artifact", show=False)
print(type(ax).__name__, ax.axison)
ax = dtcc.datasets.smoke.plot(bounds=bounds, mode="plot", show=False)
print(type(ax).__name__, ax.get_xlabel(), ax.get_ylabel())
PY
```

Expected result:

- Artifact call prints an axes type and `False` or otherwise confirms axes are off.
- Plot call prints an axes type and non-empty axis labels.

Status: completed

### Milestone 3: Implement rich preview mode

Expected changes:

- Implement `mode="preview"` as a full-figure composition.
- Return a Matplotlib `Figure`.
- Reject supplied `ax` with a clear `ValueError`.
- Use smoke presentation metadata to populate headline, summary, narrative cards, legend, annotations, key facts, and limitations.
- Draw main visual as speed slice plus streamline overlay.
- Use DTCC dark theme colors, rounded panels, integrated legend ramp, glow streamlines, and concise callouts.
- Ensure `show=False` avoids opening a window and `show=True` delegates to existing show behavior or equivalent.

Verification:

```bash
python -m pytest tests -k "smoke and preview"
```

Manual check:

```bash
python - <<'PY'
import matplotlib
matplotlib.use("Agg")
import dtcc_core as dtcc
bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
fig = dtcc.datasets.smoke.plot(bounds=bounds, mode="preview", show=False)
fig.savefig("/tmp/smoke_preview.png", dpi=120)
print(type(fig).__name__, len(fig.axes))
PY
```

Expected result:

- A non-empty `/tmp/smoke_preview.png` is created.
- The figure has multiple axes or artists consistent with a composed preview layout.
- The preview is visually polished enough for a developer/demo audience.

Status: completed

### Milestone 4: Update demo and documentation comments

Expected changes:

- Update `demos/smoke.py` to demonstrate:
  - artifact mode;
  - preview mode;
  - plot mode;
  - legacy `ax` behavior if useful.
- Keep comments clear about intended usage:
  - artifact is for table/export consumption;
  - preview is a Python preview of the table-style experience;
  - plot is ordinary technical plotting.
- Add minimal inline documentation/docstrings for `mode` and legacy `profile` behavior.
- If there is a user-facing docs location for datasets, add a small example there. Do not create broad docs work if no clear location exists.

Verification:

```bash
python demos/smoke.py
```

Expected result:

- Demo runs without errors in a local graphical environment.
- Non-graphical environments may need Agg; document if needed.

Status: completed

## Verification plan

Targeted checks:

```bash
python -m pytest tests -k "smoke"
python -m pytest tests -k "dataset"
```

Focused manual smoke test with non-interactive backend:

```bash
python - <<'PY'
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

fig_preview = dtcc.datasets.smoke.plot(bounds=bounds, show=False)
fig_preview.savefig("/tmp/smoke_preview_default.png", dpi=120)
print("preview", type(fig_preview).__name__, len(fig_preview.axes))

ax_artifact = dtcc.datasets.smoke.plot(bounds=bounds, mode="artifact", show=False)
ax_artifact.figure.savefig("/tmp/smoke_artifact.png", dpi=120)
print("artifact", type(ax_artifact).__name__, ax_artifact.axison)

ax_plot = dtcc.datasets.smoke.plot(bounds=bounds, mode="plot", show=False)
ax_plot.figure.savefig("/tmp/smoke_plot.png", dpi=120)
print("plot", type(ax_plot).__name__, ax_plot.get_xlabel(), ax_plot.get_ylabel())

fig, ax = plt.subplots()
returned = dtcc.datasets.smoke.plot(bounds=bounds, ax=ax, show=False)
print("ax identity", returned is ax)
PY
```

Expected results:

- `/tmp/smoke_preview_default.png` exists and is a rich composed preview.
- `/tmp/smoke_artifact.png` exists and contains only the visual artifact.
- `/tmp/smoke_plot.png` exists and looks like an ordinary Python plot.
- `ax identity True` is printed.

Broader checks, if practical:

```bash
python -m pytest
```

Expected results:

- All targeted tests pass.
- Full test suite passes or any unrelated failures are documented in implementation notes.

## Risks and edge cases

- Preview mode may become too smoke-specific. Keep smoke polished, but place only genuinely reusable helpers in generic plotting modules.
- Returning `Figure` for preview and `Axes` for artifact/plot is a small API inconsistency. Document it clearly. This is acceptable because preview owns a full composed layout while the other modes render into an axes.
- Existing callers without `ax` may expect an `Axes` from `datasets.smoke.plot(...)`. This behavior intentionally changes to support the new default preview experience. Existing callers using `ax` are protected by the `ax` default-to-plot rule.
- `profile` and `mode` can conflict. Do not silently choose one.
- Text wrapping in Matplotlib can be fragile. Prefer short copy and helper functions that manually wrap text to predictable line lengths.
- Preview annotation positions are relative to the visual card, not physical CRS coordinates. Label this clearly in presentation metadata and validate accordingly.
- Overlaying slice and streamlines could obscure the scalar field. Tune alpha, linewidth, and glow so both layers remain readable.
- The artifact mode must not accidentally include preview text. Tests should inspect axes/title/text artists where possible.
- Transparent backgrounds may conflict with dark preview styling. For preview mode, if `transparent=True` is requested, either support it intentionally or fail with a clear error. Do not silently produce unreadable white/missing text.
- MP4 export is currently table/artifact-oriented. Do not attempt to animate the full preview panel unless explicitly scoped later.
- The generic `DatasetDescriptor.create_context()` change could affect other datasets. Keep the merge conservative and test fallback behavior when descriptor has no presentation payload.
- If Matplotlib layout differs across versions, avoid pixel-perfect tests.

## Implementation notes

Codex should append notes here as work proceeds.

Use this section for:

- discoveries that change the plan;
- deviations from the original approach;
- decisions made during implementation;
- commands run and important results;
- risks that remain.

### Notes

- 2026-07-01: Plan created.
- 2026-07-01: Added smoke plot mode resolution for artifact, preview, and plot. Legacy profile aliases are accepted only when compatible and now fail loudly on conflicts.
- 2026-07-01: Added structured smoke presentation metadata to `SmokeDataset.describe()` and updated `DatasetDescriptor.create_context()` to consume descriptor-level presentation payloads when provided.
- 2026-07-01: Implemented smoke artifact and preview overview rendering as a speed slice with streamline overlay. Plot mode keeps ordinary Matplotlib axes behavior and supplied-axes compatibility.
- 2026-07-01: Added focused smoke plot-mode tests and updated existing smoke/context tests for the new structured presentation. Initial targeted smoke run passed: `../venv/bin/python -m pytest tests/datasets/test_smoke_plot_modes.py tests/datasets/test_smoke_dataset.py tests/datasets/test_dataset_context.py -k "smoke"`.
- 2026-07-01: Manual non-interactive render check passed and wrote `/tmp/smoke_preview_default.png`, `/tmp/smoke_artifact.png`, and `/tmp/smoke_plot.png`; preview and artifact images were visually inspected.
- 2026-07-01: Verification passed: `MPLCONFIGDIR=/tmp/dtcc-matplotlib-tests ../venv/bin/python -m pytest tests -k "smoke"` (`80 passed`) and `MPLCONFIGDIR=/tmp/dtcc-matplotlib-tests ../venv/bin/python -m pytest tests -k "dataset"` (`512 passed`).
- 2026-07-01: Demo verification passed with `MPLBACKEND=Agg MPLCONFIGDIR=/tmp/dtcc-matplotlib-tests ../venv/bin/python demos/smoke.py`; Agg emitted the expected non-interactive `plt.show()` warning.

## Decision log

Record important implementation decisions.

| Date | Decision | Reason |
|---|---|---|
| 2026-07-01 | Use one public `mode` option with values `artifact`, `preview`, and `plot`. | Avoids confusing `profile + view` combinations and maps directly to real use cases. |
| 2026-07-01 | Keep `profile` temporarily as a legacy alias. | Preserves existing smoke export/plot behavior while introducing the simpler API. |
| 2026-07-01 | Default `plot()` without `ax` to `preview`; default `plot(..., ax=ax)` to `plot`. | Provides the desired rich dataset preview while preserving common Matplotlib subplot usage. |
| 2026-07-01 | Preview mode returns `Figure`; artifact/plot modes return `Axes`. | Preview owns a full figure layout, while artifact/plot can render into a single axes. |
| 2026-07-01 | Artifact mode contains no narrative, axes, title, or frontend-owned UX. | Tangible-table frontend should own final layout, legends, and narrative placement. |
| 2026-07-01 | Smoke preview uses a speed-slice background plus streamline overlay. | This is the strongest single visual for demonstrating field magnitude and flow direction. |

## Final review checklist

Before this task is accepted:

- [x] Acceptance criteria are satisfied.
- [x] Required data/configuration fails loudly when missing or invalid.
- [x] No silent fallbacks or placeholder defaults were introduced.
- [x] Human-facing CLI behavior is simple for the common case, if applicable.
- [x] Tests were added or updated for changed behavior.
- [x] Verification commands were run, or limitations were documented.
- [x] No unrelated refactors or broad rewrites were introduced.
- [x] Public APIs remain compatible unless the plan explicitly changes them.
- [x] Security, authorization, data integrity, and migration risks were considered.
- [x] No known blocking issues remain.

## Done condition

The task is done when the acceptance criteria are met, relevant verification has passed or limitations are documented, and review finds no blocking correctness, safety, test, or maintainability issues.
