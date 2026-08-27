"""Small, shared BigWig operations used by the CUT&RUN/CUT&Tag workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

try:
    import pyBigWig
except ImportError:  # pragma: no cover - dependency validation reports this first
    pyBigWig = None


def _require_pybigwig() -> None:
    if pyBigWig is None:
        raise RuntimeError("pyBigWig is required for CUT&RUN/CUT&Tag BigWig comparisons")


def _finite_values(handle, chrom: str, start: int, end: int) -> np.ndarray:
    if end <= start:
        return np.empty(0, dtype=np.float64)
    values = np.asarray(handle.values(chrom, start, end, numpy=True), dtype=np.float64)
    return values[np.isfinite(values)]


def interval_max(handle, chrom: str, start: int, end: int) -> float:
    """Return the maximum finite value in an interval, or zero when it is empty."""

    values = _finite_values(handle, chrom, start, end)
    return float(np.max(values)) if values.size else 0.0


def interval_positive_area(handle, chrom: str, start: int, end: int) -> float:
    """Return the base-wise area above zero in an interval."""

    values = _finite_values(handle, chrom, start, end)
    return float(np.sum(np.maximum(values, 0.0))) if values.size else 0.0


def bigwig_chroms(path: str | Path) -> dict[str, int]:
    _require_pybigwig()
    handle = pyBigWig.open(str(path))
    if handle is None:
        raise ValueError(f"Could not open BigWig: {path}")
    try:
        return {str(name): int(length) for name, length in handle.chroms().items()}
    finally:
        handle.close()


def average_bigwigs(
    paths: Sequence[str | Path],
    output_path: str | Path,
    *,
    chunk_bp: int = 1_000_000,
) -> Path:
    """Write the arithmetic base-wise mean of compatible BigWigs."""

    _require_pybigwig()
    if not paths:
        raise ValueError("At least one BigWig is required")
    if chunk_bp < 1:
        raise ValueError("chunk_bp must be positive")
    inputs = [Path(path).resolve() for path in paths]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)

    handles = [pyBigWig.open(str(path)) for path in inputs]
    if any(handle is None for handle in handles):
        for handle in handles:
            if handle is not None:
                handle.close()
        raise ValueError("Could not open one or more BigWig inputs")
    try:
        chroms = {str(name): int(length) for name, length in handles[0].chroms().items()}
        for path, handle in zip(inputs[1:], handles[1:]):
            observed = {
                str(name): int(length) for name, length in handle.chroms().items()
            }
            if observed != chroms:
                raise ValueError(
                    f"BigWig chromosome names or lengths differ: {inputs[0]} and {path}"
                )

        output = Path(output_path).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        writer = pyBigWig.open(str(output), "w")
        if writer is None:
            raise ValueError(f"Could not create BigWig: {output}")
        try:
            writer.addHeader(list(chroms.items()))
            divisor = float(len(handles))
            for chrom, length in chroms.items():
                for start in range(0, length, chunk_bp):
                    end = min(length, start + chunk_bp)
                    total = np.zeros(end - start, dtype=np.float64)
                    for handle in handles:
                        values = np.asarray(
                            handle.values(chrom, start, end, numpy=True),
                            dtype=np.float64,
                        )
                        total += np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
                    mean = total / divisor
                    writer.addEntries(
                        chrom,
                        start,
                        values=mean.astype(float).tolist(),
                        span=1,
                        step=1,
                    )
        finally:
            writer.close()
        return output
    finally:
        for handle in handles:
            handle.close()


def open_bigwigs(paths: Iterable[str | Path]) -> list[object]:
    """Open BigWigs and close already-opened handles if one input fails."""

    _require_pybigwig()
    handles: list[object] = []
    try:
        for path in paths:
            handle = pyBigWig.open(str(path))
            if handle is None:
                raise ValueError(f"Could not open BigWig: {path}")
            handles.append(handle)
        return handles
    except Exception:
        for handle in handles:
            handle.close()
        raise
