"""Durable per-invocation logs for executing NucleoSuite commands."""

from __future__ import annotations

import json
import os
import shlex
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Mapping, Sequence, TextIO


def _option_value(argv: Sequence[str], names: Sequence[str]) -> str | None:
    value: str | None = None
    for index, token in enumerate(argv):
        for name in names:
            if token == name and index + 1 < len(argv):
                value = argv[index + 1]
            elif token.startswith(name + "="):
                value = token.split("=", 1)[1]
    return value


def default_log_path(command: str, argv: Sequence[str]) -> Path:
    """Choose a deterministic output-adjacent directory for an invocation log."""
    directory_value = _option_value(argv, ("--outdir", "--output-dir"))
    prefix_value = _option_value(argv, ("--output-prefix", "--out-prefix", "-o"))
    if directory_value:
        root = Path(directory_value)
    elif prefix_value:
        root = Path(prefix_value).parent
    else:
        root = Path.cwd() / "nucleosuite_logs"
    # Microseconds keep rapid programmatic invocations in distinct files.  The
    # PID alone is not sufficient when a long-lived Python process calls the
    # CLI more than once within the same second.
    timestamp = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S.%f%z")
    safe_command = command.replace("/", "_").replace(" ", "_")
    return root / "logs" / "commands" / f"{timestamp}_{safe_command}_{os.getpid()}.log"


def serializable_parameters(parameters: Mapping[str, object] | None) -> dict[str, object]:
    output: dict[str, object] = {}
    for key, value in sorted((parameters or {}).items()):
        if key.startswith("_") or callable(value):
            continue
        if isinstance(value, Path):
            output[key] = str(value)
        elif isinstance(value, (str, int, float, bool)) or value is None:
            output[key] = value
        elif isinstance(value, (list, tuple)):
            output[key] = [str(item) if isinstance(item, Path) else item for item in value]
        elif isinstance(value, dict):
            output[key] = {str(k): str(v) for k, v in value.items()}
        else:
            output[key] = str(value)
    return output


class _Tee:
    def __init__(self, terminal: TextIO, log: TextIO) -> None:
        self.terminal = terminal
        self.log = log

    def write(self, text: str) -> int:
        self.terminal.write(text)
        self.log.write(text)
        return len(text)

    def flush(self) -> None:
        self.terminal.flush()
        self.log.flush()

    def isatty(self) -> bool:
        return bool(getattr(self.terminal, "isatty", lambda: False)())

    @property
    def encoding(self) -> str | None:
        return getattr(self.terminal, "encoding", None)


class CommandLog:
    """Mirror Python console output and record command metadata durably."""

    def __init__(
        self,
        command: str,
        argv: Sequence[str],
        *,
        version: str,
        parameters: Mapping[str, object] | None = None,
    ) -> None:
        self.command = command
        self.argv = list(argv)
        self.version = version
        self.parameters = serializable_parameters(parameters)
        self.path = default_log_path(command, argv)
        self.started = time.monotonic()
        self._handle: TextIO | None = None
        self._stdout: TextIO | None = None
        self._stderr: TextIO | None = None

    def __enter__(self) -> "CommandLog":
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._handle = self.path.open("a", encoding="utf-8", buffering=1)
        except OSError:
            fallback = Path.cwd() / "nucleosuite_logs" / self.path.name
            fallback.parent.mkdir(parents=True, exist_ok=True)
            self.path = fallback
            self._handle = self.path.open("a", encoding="utf-8", buffering=1)
        self._handle.write(
            f"[START] {datetime.now().astimezone().isoformat(timespec='seconds')}\n"
        )
        self._handle.write(f"[VERSION] {self.version}\n")
        self._handle.write(
            "[COMMAND] " + shlex.join(["nucleosuite", *self.argv]) + "\n"
        )
        self._handle.write(
            "[PARAMETERS] "
            + json.dumps(self.parameters, sort_keys=True, separators=(",", ":"))
            + "\n"
        )
        self._handle.write(f"[WORKING_DIRECTORY] {Path.cwd()}\n")
        self._stdout, self._stderr = sys.stdout, sys.stderr
        sys.stdout = _Tee(self._stdout, self._handle)  # type: ignore[assignment]
        sys.stderr = _Tee(self._stderr, self._handle)  # type: ignore[assignment]
        print(f"[{self.command}] Command log: {self.path}", flush=True)
        return self

    def finish(self, exit_code: int) -> None:
        if self._handle is None:
            return
        self._handle.write(
            f"[END] {datetime.now().astimezone().isoformat(timespec='seconds')} "
            f"exit_code={int(exit_code)} elapsed_seconds={time.monotonic() - self.started:.3f}\n"
        )
        self._handle.flush()

    def __exit__(self, exc_type, exc, traceback) -> None:
        if exc_type is not None:
            self.finish(1)
        if self._stdout is not None:
            sys.stdout = self._stdout
        if self._stderr is not None:
            sys.stderr = self._stderr
        if self._handle is not None:
            self._handle.close()
