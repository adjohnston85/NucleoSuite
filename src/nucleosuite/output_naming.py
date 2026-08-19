"""Stable, collision-resistant output-prefix helpers.

Analysis parameters belong in automatic filenames because changing them changes
the data represented by every related table and plot.  Presentation-only
settings (DPI, colours, fonts, and tick spacing) deliberately do not.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable


def compact_parameter(value: object) -> str:
    """Return a concise filename-safe representation of an option value."""

    if value is None:
        return "none"
    if isinstance(value, bool):
        return "on" if value else "off"
    numeric = isinstance(value, (int, float)) and not isinstance(value, bool)
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if math.isinf(value):
            return "inf" if value > 0 else "neginf"
        text = str(int(value)) if value.is_integer() else f"{value:.12g}"
    else:
        text = str(value)
    text = text.strip().lower()
    if numeric and text.startswith("-"):
        text = "neg" + text[1:]
    text = text.replace("+", "")
    text = text.replace(".", "p").replace("%", "pct")
    text = re.sub(r"[^a-z0-9]+", "-", text).strip("-")
    return text or "none"


def parameterized_prefix(
    base: str | Path,
    parameters: Iterable[tuple[str, object]],
    *,
    marker: str | None = None,
) -> Path:
    """Append ordered analysis parameter tokens to *base* exactly once."""

    path = Path(base)
    tokens = [f"{key}{compact_parameter(value)}" for key, value in parameters]
    if not tokens:
        return path
    sentinel = marker or tokens[0]
    if re.search(rf"(?:^|_){re.escape(sentinel)}(?:_|$)", path.name):
        return path
    return path.with_name(path.name + "_" + "_".join(tokens))


def parameter_range(start: object, end: object) -> str:
    """Format a signed interval as one compact parameter value."""

    return f"{compact_parameter(start)}to{compact_parameter(end)}"
