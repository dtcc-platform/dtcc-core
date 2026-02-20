# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

"""Custom Rich logging handler for DTCC."""

from datetime import datetime
from logging import LogRecord

from rich.console import Console
from rich.logging import RichHandler
from rich.text import Text
from rich.traceback import Traceback
from rich.highlighter import NullHighlighter

# Maps log level -> (level_style, message_style)
LEVEL_STYLES = {
    "DEBUG": ("dim", "dim"),
    "INFO": ("bright_blue", ""),  # empty = default
    "WARNING": ("yellow", "yellow"),
    "ERROR": ("red bold", "red"),
    "CRITICAL": ("white on red", "red bold"),
}


class LoggingHandler(RichHandler):
    """Custom logging handler with DTCC-style formatting.

    Output format:
        [13:24:20] [dtcc-core] [info] Message text here
    """

    def __init__(
        self,
        source_name: str = "dtcc-core",
        console: Console = None,
        **kwargs,
    ):
        # Disable RichHandler's built-in formatting
        super().__init__(
            console=console,
            show_time=False,
            show_level=False,
            show_path=False,
            rich_tracebacks=True,
            markup=False,
            highlighter=NullHighlighter(),
            **kwargs,
        )
        self.source_name = source_name

    def emit(self, record: LogRecord) -> None:
        """Emit a log record with DTCC formatting."""
        try:
            # Get styles for this level
            level_style, msg_style = LEVEL_STYLES.get(
                record.levelname, ("", "")
            )

            # Format timestamp
            time_str = datetime.fromtimestamp(record.created).strftime("%H:%M:%S")

            # Format level name (lowercase, no padding)
            level_name = record.levelname.lower()

            # Get the formatted message
            message = record.getMessage()

            # Handle multiline messages
            lines = message.split("\n")

            for line_text in lines:
                line = Text.assemble(
                    ("[", "dim"),
                    (time_str, "dim"),
                    ("] ", "dim"),
                    ("[", "dim"),
                    (self.source_name, "dim"),
                    ("] ", "dim"),
                    ("[", "dim"),
                    (level_name, level_style),
                    ("] ", "dim"),
                    (line_text, msg_style or ""),
                )
                self.console.print(line)

            # Handle exceptions/tracebacks
            if record.exc_info:
                tb = Traceback.from_exception(
                    record.exc_info[0],
                    record.exc_info[1],
                    record.exc_info[2],
                    show_locals=self.tracebacks_show_locals,
                )
                self.console.print(tb)

        except Exception:
            self.handleError(record)
