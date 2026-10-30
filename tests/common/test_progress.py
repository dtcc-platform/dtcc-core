"""Tests for progress tracking monotonic guarantee."""

import os
import sys

# Import directly from the source tree to pick up local changes,
# bypassing the dtcc_core __init__.py which loads native extensions.
sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "..", "dtcc_core", "common")
)
from progress import ProgressTracker, report_progress


class TestMonotonicProgress:
    """Progress within a phase must never decrease."""

    def test_lower_percent_does_not_decrease_progress(self):
        with ProgressTracker(phases={"work": 1.0}, mode="silent") as tracker:
            with tracker.phase("work"):
                report_progress(percent=80)
                phase = tracker.state.phases["work"]
                assert phase.progress == 0.8

                report_progress(percent=60)
                assert phase.progress == 0.8  # must not go backward

                report_progress(percent=90)
                assert phase.progress == 0.9  # higher value still accepted

    def test_nested_call_pattern_no_regression(self):
        """Simulate: parent reports 10%, child reports 0%->100%, parent reports 30%."""
        with ProgressTracker(phases={"build": 1.0}, mode="silent") as tracker:
            with tracker.phase("build"):
                # Parent starts
                report_progress(percent=10)
                phase = tracker.state.phases["build"]
                assert phase.progress == 0.1

                # Child drives to 100%
                report_progress(percent=0)
                assert phase.progress == 0.1  # must not drop
                report_progress(percent=50)
                assert phase.progress == 0.5
                report_progress(percent=100)
                assert phase.progress == 1.0

                # Parent resumes at 30% — must not regress
                report_progress(percent=30)
                assert phase.progress == 1.0

    def test_overall_percent_monotonic_across_phases(self):
        with ProgressTracker(phases={"a": 0.5, "b": 0.5}, mode="silent") as tracker:
            prev = 0.0
            with tracker.phase("a"):
                for pct in [0, 25, 50, 75, 100]:
                    report_progress(percent=pct)
                    current = tracker.state.overall_percent
                    assert current >= prev, (
                        f"overall_percent went backward: {prev} -> {current}"
                    )
                    prev = current

            with tracker.phase("b"):
                for pct in [0, 25, 50, 75, 100]:
                    report_progress(percent=pct)
                    current = tracker.state.overall_percent
                    assert current >= prev, (
                        f"overall_percent went backward: {prev} -> {current}"
                    )
                    prev = current


def test_nested_progress_ranges_restore_parent_and_bound_child():
    import pytest
    from progress import progress_range

    with ProgressTracker(phases={"mesh": 1.0}, mode="silent") as tracker:
        with tracker.phase("mesh"):
            with progress_range(20, 80):
                report_progress(percent=50)
                assert tracker.state.overall_percent == pytest.approx(50)
                with progress_range(50, 75):
                    report_progress(percent=50)
                    assert tracker.state.overall_percent == pytest.approx(57.5)
                    report_progress(increment=1000)
                    assert tracker.state.overall_percent == pytest.approx(65)
                report_progress(percent=90)
                assert tracker.state.overall_percent == pytest.approx(74)
            assert tracker.state.overall_percent == pytest.approx(80)
            report_progress(percent=95)
            assert tracker.state.overall_percent == pytest.approx(95)


def test_progress_range_failure_does_not_complete_and_new_tracker_is_independent():
    import pytest
    from progress import progress_range

    with ProgressTracker(phases={"mesh": 1.0}, mode="silent") as outer:
        with outer.phase("mesh"):
            with pytest.raises(RuntimeError):
                with progress_range(20, 80):
                    report_progress(percent=25)
                    with ProgressTracker(phases={"other": 1.0}, mode="silent") as inner:
                        with inner.phase("other"):
                            report_progress(percent=50)
                            assert inner.state.overall_percent == 50
                    raise RuntimeError("controlled failure")
            assert outer.state.overall_percent == pytest.approx(35)
            report_progress(percent=40)
            assert outer.state.overall_percent == pytest.approx(40)
