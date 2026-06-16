from pathlib import Path

import dtcc_core as dtcc


OUTPUT_ROOT = Path(__file__).resolve().parent / "output"


# Bounds are in EPSG:3006 / SWEREF 99 TM.
AREAS = (
    {
        "name": "verkstadsgatan_3",
        "bounds": dtcc.Bounds(
            xmin=324861,
            ymin=6374485,
            xmax=325668,
            ymax=6375155,
        ),
    },
    {
        "name": "turkosvagen_4",
        "bounds": dtcc.Bounds(
            xmin=323963,
            ymin=6374269,
            xmax=324426,
            ymax=6374687,
        ),
    },
    {
        "name": "alvsakersvagen_500",
        "bounds": dtcc.Bounds(
            xmin=330476,
            ymin=6381158,
            xmax=331668,
            ymax=6382045,
        ),
    },
)
