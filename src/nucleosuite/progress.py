"""Small, rate-limited progress reporter shared by long-running commands."""

from __future__ import annotations

import time
import sys
from dataclasses import dataclass, field
from typing import TextIO


@dataclass
class ProgressReporter:
    command: str
    quiet: bool = False
    minimum_interval_seconds: float = 5.0
    stream: TextIO = field(default_factory=lambda: sys.stderr)
    _started: float = field(default_factory=time.monotonic)
    _last_update: float = field(default=0.0)

    def emit(self, message: str, *, force: bool = True) -> None:
        now = time.monotonic()
        if self.quiet:
            return
        if not force and now - self._last_update < self.minimum_interval_seconds:
            return
        self._last_update = now
        print(f"[{self.command}] {message}", file=self.stream, flush=True)

    def stage(self, message: str) -> None:
        self.emit(message)

    def file_start(self, label: str, path: object) -> None:
        self.emit(f"Reading {label}: {path}")

    def file_progress(self, label: str, records: int, contig: str) -> None:
        """Emit an optional rate-limited record update.

        Commands should normally prefer one message per new contig.  This
        method remains available for operations that cannot expose contig
        transitions, but is deliberately rate limited.
        """
        self.emit(
            f"{label}: {contig}; {records:,} records read",
            force=False,
        )

    def reading_contig(self, label: str, contig: str) -> None:
        """Report the current input contig once, without cumulative counters."""
        self.emit(f"Reading {label} contig: {contig}")

    def file_complete(self, label: str, path: object, records: int) -> None:
        self.emit(f"Completed {label}: {path}; {records:,} records")

    def contig(
        self,
        stage: str,
        contig: str,
        index: int,
        total: int,
        records: int | None = None,
    ) -> None:
        suffix = "" if records is None else f"; {records:,} query records"
        self.emit(f"{stage}: {contig} ({index}/{total}){suffix}")

    def complete(self, message: str = "Completed") -> None:
        self.emit(f"{message}; elapsed {time.monotonic() - self._started:.1f} s")
