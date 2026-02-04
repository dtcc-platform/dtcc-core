# dtcc_core/progress.py
"""
Progress tracking system for DTCC.

Usage at Dataset level:

    class MyDataset(DatasetDescriptor):
        def build(self, args):
            with ProgressTracker(phases={"download": 0.1, "process": 0.3, "mesh": 0.6}) as progress:
                with progress.phase("download"):
                    data = download_data()

                with progress.phase("process"):
                    processed = process_data(data)

                with progress.phase("mesh"):
                    result = build_mesh(processed)

            return result

Usage in nested functions (no changes needed to function signatures):

    def some_inner_function():
        # This "just works" if called within a ProgressTracker context
        report_progress(0.5, message="Halfway through inner function")

"""

import sys
import time
import json
import threading
import functools
import contextvars
from typing import Optional, Dict, Iterator, Callable, Any, List
from dataclasses import dataclass, field
from contextlib import contextmanager
from enum import Enum, auto

from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
    TimeRemainingColumn,
    TaskID,
)
from rich.console import Console

# Shared console for progress and logging integration
# Using stderr to keep stdout clean for actual program output
_console = Console(stderr=True)


def get_console() -> Console:
    """Get the shared rich Console instance for logging integration."""
    return _console


class ProgressMode(Enum):
    TERMINAL = auto()  # Human-readable terminal output
    JSON = auto()       # Machine-readable JSON lines
    SILENT = auto()     # No output (for nested trackers)
    CALLBACK = auto()   # Custom callback function


# ============================================================
# Context Variable for Propagation
# ============================================================

_current_progress: contextvars.ContextVar[Optional['ProgressTracker']] = \
    contextvars.ContextVar('dtcc_progress', default=None)


def get_progress() -> Optional['ProgressTracker']:
    """Get the current progress tracker from context."""
    return _current_progress.get()


def report_progress(
    percent: float = None,
    message: str = None,
    increment: float = None,
    current: int = None,
    total: int = None,
):
    """
    Report progress from anywhere in the call stack.
    Safe to call even if no tracker is active (becomes a no-op).

    Args:
        percent: Progress within current phase (0-100)
        message: Status message to display
        increment: Increment progress by this percentage
        current: Current item number (alternative to percent)
        total: Total items (used with current)

    Examples:
        report_progress(percent=50, message="Halfway done")
        report_progress(current=5, total=10)  # 50%
        report_progress(increment=10)  # Add 10%
    """
    tracker = _current_progress.get()
    if tracker is None:
        return

    # Convert current/total to percent
    if current is not None and total is not None and total > 0:
        percent = (current / total) * 100

    tracker._update_phase_progress(percent=percent, increment=increment, message=message)


# ============================================================
# Progress State
# ============================================================

@dataclass
class PhaseInfo:
    """Information about a single phase."""
    name: str
    weight: float
    progress: float = 0.0  # 0-1 within this phase
    message: str = ""
    started: bool = False
    completed: bool = False


@dataclass
class ProgressState:
    """Complete state of progress tracking."""
    phases: Dict[str, PhaseInfo] = field(default_factory=dict)
    current_phase: Optional[str] = None
    start_time: float = field(default_factory=time.time)
    message: str = ""

    # For simple (non-phased) progress
    current: int = 0
    total: int = 0

    @property
    def overall_percent(self) -> float:
        """Calculate overall progress from weighted phases."""
        if not self.phases:
            # Simple mode
            if self.total <= 0:
                return 0.0
            return (self.current / self.total) * 100

        total = 0.0
        for phase in self.phases.values():
            if phase.completed:
                total += phase.weight * 100
            elif phase.started:
                total += phase.weight * phase.progress * 100
        return total

    @property
    def eta_seconds(self) -> Optional[float]:
        """Estimate time remaining."""
        elapsed = time.time() - self.start_time
        pct = self.overall_percent
        if pct <= 0 or elapsed <= 0:
            return None
        return elapsed * (100 - pct) / pct


# ============================================================
# Main Progress Tracker
# ============================================================

class ProgressTracker:
    """
    Thread-safe progress tracker with phase support and context propagation.

    Args:
        phases: Dict mapping phase names to their weights (should sum to 1.0)
        mode: Output mode (auto-detected if not specified)
        output: Output stream (default: stderr)
        callback: Custom callback function for CALLBACK mode
        min_update_interval: Minimum seconds between renders (prevents flickering)

    Examples:
        # Phased progress
        with ProgressTracker(phases={"load": 0.2, "process": 0.8}) as p:
            with p.phase("load"):
                load_data()
            with p.phase("process"):
                process_data()

        # Simple progress
        with ProgressTracker(total=100) as p:
            for i in range(100):
                p.update(current=i)
    """

    def __init__(
        self,
        phases: Dict[str, float] = None,
        total: float = 0.0,
        mode: str = "auto",
        output = None,
        callback: Callable[[dict], None] = None,
        min_update_interval: float = 0.05,
    ):
        self.state = ProgressState(total=total)
        self.output = output or sys.stderr
        self.callback = callback
        self.min_update_interval = min_update_interval
        self._lock = threading.Lock()
        self._token = None
        self._last_render_time = 0

        # Rich progress components
        self._rich_progress: Optional[Progress] = None
        self._rich_task_id: Optional[TaskID] = None  # Single unified task
        self._console = _console

        # Set up phases
        if phases:
            # Normalize weights to sum to 1.0
            total_weight = sum(phases.values())
            self.state.phases = {
                name: PhaseInfo(name=name, weight=w/total_weight)
                for name, w in phases.items()
            }

        # Determine output mode
        if mode == "auto":
            if callback:
                self._mode = ProgressMode.CALLBACK
            else:
                try:
                    self._mode = ProgressMode.TERMINAL if self.output.isatty() else ProgressMode.JSON
                except AttributeError:
                    self._mode = ProgressMode.JSON
        else:
            self._mode = {
                "terminal": ProgressMode.TERMINAL,
                "json": ProgressMode.JSON,
                "silent": ProgressMode.SILENT,
                "callback": ProgressMode.CALLBACK,
            }.get(mode, ProgressMode.JSON)

    # ---- Context Manager ----

    def __enter__(self):
        self._token = _current_progress.set(self)
        self.state.start_time = time.time()

        # Start rich Progress for terminal mode
        if self._mode == ProgressMode.TERMINAL:
            self._rich_progress = Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                TimeRemainingColumn(),
                console=self._console,
                transient=True,
            )
            self._rich_progress.start()

            # Create a single task for overall progress
            self._rich_task_id = self._rich_progress.add_task(
                "Starting...", total=100, visible=True
            )

        return self

    def __exit__(self, *args):
        if self._token:
            _current_progress.reset(self._token)
        self._finish()

    # ---- Phase Management ----

    @contextmanager
    def phase(self, name: str, message: str = None):
        """
        Context manager for entering a progress phase.

        Args:
            name: Phase name (must match a key in phases dict)
            message: Optional status message
        """
        if name not in self.state.phases:
            # Allow ad-hoc phases (equal weight)
            self.state.phases[name] = PhaseInfo(name=name, weight=1.0)

        phase_info = self.state.phases[name]
        phase_info.started = True
        phase_info.message = message or name
        self.state.current_phase = name
        self.state.message = message or name

        # Update the single progress bar with current phase
        if self._rich_progress and self._rich_task_id is not None:
            self._rich_progress.update(
                self._rich_task_id,
                description=message or name,
                completed=self.state.overall_percent,
                refresh=True,
            )

        self._render()

        try:
            yield
        finally:
            phase_info.progress = 1.0
            phase_info.completed = True

            # Update overall progress on phase completion
            if self._rich_progress and self._rich_task_id is not None:
                self._rich_progress.update(
                    self._rich_task_id,
                    completed=self.state.overall_percent,
                    refresh=True,
                )

            self._render()

    def _update_phase_progress(
        self,
        percent: float = None,
        increment: float = None,
        message: str = None
    ):
        """Update progress within current phase."""
        with self._lock:
            if self.state.current_phase and self.state.current_phase in self.state.phases:
                phase = self.state.phases[self.state.current_phase]

                if increment is not None:
                    phase.progress = min(1.0, phase.progress + increment / 100)
                elif percent is not None:
                    phase.progress = min(1.0, percent / 100)

                if message:
                    phase.message = message
                    self.state.message = message

                # Update the single progress bar with overall progress
                if self._rich_progress and self._rich_task_id is not None:
                    update_kwargs = {"completed": self.state.overall_percent, "refresh": True}
                    if message:
                        update_kwargs["description"] = message
                    self._rich_progress.update(self._rich_task_id, **update_kwargs)

            elif message:
                self.state.message = message

            self._render()

    # ---- Simple (Non-Phased) Updates ----

    def update(
        self,
        current: int = None,
        total: int = None,
        increment: int = None,
        message: str = None
    ):
        """Update simple progress (non-phased mode)."""
        with self._lock:
            if total is not None:
                self.state.total = total
            if increment is not None:
                self.state.current += increment
            elif current is not None:
                self.state.current = current
            if message:
                self.state.message = message

            # Update rich progress for simple mode
            if self._rich_progress and self._rich_task_id is not None:
                update_kwargs = {"completed": self.state.overall_percent, "refresh": True}
                if message:
                    update_kwargs["description"] = message
                self._rich_progress.update(self._rich_task_id, **update_kwargs)

            self._render()

    def set_total(self, total: int):
        """Set total count for simple progress."""
        with self._lock:
            self.state.total = total

    # ---- Iterable Tracking ----

    def track(self, iterable, total: int = None, message: str = None) -> Iterator:
        """
        Wrap an iterable to track progress automatically.

        Args:
            iterable: Any iterable
            total: Total count (inferred from len() if not provided)
            message: Status message template (can include {current}, {total})
        """
        if total is None:
            try:
                total = len(iterable)
            except TypeError:
                iterable = list(iterable)
                total = len(iterable)

        self.state.total = total

        for i, item in enumerate(iterable):
            msg = message.format(current=i+1, total=total) if message and '{' in message else message
            self.update(current=i, message=msg)
            yield item

        self.update(current=total)

    # ---- C++ Integration ----

    def make_callback(self, phase_name: str = None) -> Callable[[int, int, str], None]:
        """
        Create a callback function to pass to C++ code.

        Args:
            phase_name: If provided, updates are scoped to this phase

        Returns:
            Callable with signature (current, total, message) -> None
        """
        def callback(current: int, total: int, message: str = ""):
            if phase_name and phase_name in self.state.phases:
                percent = (current / total * 100) if total > 0 else 0
                self._update_phase_progress(percent=percent, message=message)
            else:
                self.update(current=current, total=total, message=message)

        return callback

    # ---- Rendering ----

    def _render(self):
        """Render progress to output."""
        # Rate limiting
        now = time.time()
        if now - self._last_render_time < self.min_update_interval:
            return
        self._last_render_time = now

        # Rich handles terminal rendering automatically
        if self._mode == ProgressMode.TERMINAL:
            # Rich Progress handles this - no additional rendering needed
            pass
        elif self._mode == ProgressMode.JSON:
            self._render_json()
        elif self._mode == ProgressMode.CALLBACK and self.callback:
            self.callback(self._get_state_dict())

    def _render_json(self):
        """Render machine-readable JSON output."""
        data = self._get_state_dict()
        # Use magic prefix for easy grep
        print(f"##PROGRESS##{json.dumps(data)}##", file=self.output, flush=True)

    def _get_state_dict(self) -> dict:
        """Get current state as dictionary."""
        return {
            "type": "progress",
            "percent": round(self.state.overall_percent, 2),
            "eta_seconds": round(self.state.eta_seconds) if self.state.eta_seconds else None,
            "eta_formatted": self._format_eta(self.state.eta_seconds),
            "message": self.state.message,
            "phase": self.state.current_phase,
            "phases": {
                name: {
                    "weight": round(p.weight, 3),
                    "progress": round(p.progress * 100, 1),
                    "completed": p.completed,
                    "message": p.message,
                }
                for name, p in self.state.phases.items()
            },
            "elapsed": round(time.time() - self.state.start_time, 2),
        }

    def _finish(self):
        """Clean up after progress completes."""
        if self._mode == ProgressMode.TERMINAL:
            # Stop rich progress
            if self._rich_progress:
                self._rich_progress.stop()
                self._rich_progress = None
        elif self._mode == ProgressMode.JSON:
            # Emit completion event
            data = self._get_state_dict()
            data["type"] = "progress_complete"
            data["percent"] = 100.0
            print(f"##PROGRESS##{json.dumps(data)}##", file=self.output, flush=True)

    @staticmethod
    def _format_eta(seconds: Optional[float]) -> str:
        """Format ETA as human-readable string."""
        if seconds is None:
            return "calculating..."
        if seconds < 0:
            return "almost done"
        if seconds < 60:
            return f"{int(seconds)}s"
        if seconds < 3600:
            return f"{int(seconds // 60)}m {int(seconds % 60)}s"
        if seconds < 86400:
            return f"{int(seconds // 3600)}h {int((seconds % 3600) // 60)}m"
        days = int(seconds // 86400)
        hours = int((seconds % 86400) // 3600)
        if days > 7:
            weeks = days // 7
            days = days % 7
            return f"{weeks}w {days}d"
        return f"{days}d {hours}h"


# ============================================================
# Decorator for Functions
# ============================================================

def with_progress(
    phases: Dict[str, float] = None,
    total: int = None,
    message: str = None,
):
    """
    Decorator that sets up progress tracking for a function.

    If a parent tracker exists, creates a sub-phase instead of new tracker.

    Args:
        phases: Phase weights for this function
        total: Total items for simple progress
        message: Default status message
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            parent = get_progress()

            if parent is None:
                # No parent - create new tracker
                with ProgressTracker(phases=phases, total=total) as tracker:
                    if message:
                        tracker.state.message = message
                    return func(*args, **kwargs)
            else:
                # Has parent - just run (progress reports go to parent)
                return func(*args, **kwargs)

        return wrapper
    return decorator
