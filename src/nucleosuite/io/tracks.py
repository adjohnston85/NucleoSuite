"""BigWig and compressed-WIG writers shared by all track commands."""

from __future__ import annotations

import gzip
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence

import numpy as np

try:
    import pyBigWig
except ImportError:  # pragma: no cover - handled with a clear runtime error
    pyBigWig = None

from nucleosuite.scoring.basic_tracks import SPARSE_TRACKS

OUTPUT_FORMATS = ("bigwig", "wiggz", "both", "none")


class TrackHandleMap(dict):
    """Handle mapping with atomic partial/final output bookkeeping."""

    def __init__(self) -> None:
        super().__init__()
        self.paths: dict[str, tuple[Path, Path]] = {}


def require_pybigwig() -> None:
    if pyBigWig is None:
        raise RuntimeError(
            "BigWig output requires pyBigWig. Install it with "
            "'conda install -c bioconda pybigwig' or 'pip install pyBigWig'."
        )


def build_bigwig_header(
    references: Sequence[str],
    lengths: Sequence[int],
    selected_contigs: Iterable[str],
):
    selected = set(selected_contigs)
    return [
        (name, int(length))
        for name, length in zip(references, lengths)
        if name in selected
    ]


def open_track_handles(
    output_prefix: str,
    tracks: Sequence[str],
    output_format: str,
    bigwig_header,
):
    if output_format not in OUTPUT_FORMATS:
        raise ValueError(f"Unknown output format: {output_format}")
    bigwig_handles = TrackHandleMap()
    wig_handles = TrackHandleMap()
    if output_format in {"bigwig", "both"}:
        require_pybigwig()
        for track in dict.fromkeys(tracks):
            final_path = Path(f"{output_prefix}_{track}.bw")
            partial_path = final_path.with_name(final_path.name + ".partial")
            partial_path.unlink(missing_ok=True)
            final_path.with_name(final_path.name + ".complete.json").unlink(
                missing_ok=True
            )
            handle = pyBigWig.open(str(partial_path), "w")
            handle.addHeader(bigwig_header)
            bigwig_handles[track] = handle
            bigwig_handles.paths[track] = (partial_path, final_path)
    if output_format in {"wiggz", "both"}:
        for track in dict.fromkeys(tracks):
            final_path = Path(f"{output_prefix}_{track}.wig.gz")
            partial_path = final_path.with_name(final_path.name + ".partial")
            partial_path.unlink(missing_ok=True)
            final_path.with_name(final_path.name + ".complete.json").unlink(
                missing_ok=True
            )
            wig_handles[track] = gzip.open(partial_path, "wt")
            wig_handles.paths[track] = (partial_path, final_path)
    return bigwig_handles, wig_handles


def remove_stale_track_outputs(output_prefix: str, tracks: Iterable[str]) -> None:
    for track in tracks:
        for suffix in ("bw", "wig.gz", "bedGraph"):
            path = f"{output_prefix}_{track}.{suffix}"
            partial = path + ".partial"
            if os.path.exists(partial):
                os.remove(partial)
            marker = path + ".complete.json"
            if os.path.exists(marker):
                os.remove(marker)


def write_bigwig_tracks(
    scores,
    contig: str,
    adjusted_start: int,
    original_start: int,
    original_end: int,
    handles: Mapping,
    tracks: Sequence[str],
) -> None:
    for track in tracks:
        if track not in scores or track not in handles:
            continue
        array = np.asarray(scores[track][0][2], dtype=np.float64)
        left = original_start - adjusted_start
        right = original_end - adjusted_start
        values = array[left:right]
        if values.size == 0:
            continue

        handle = handles[track]
        if track in SPARSE_TRACKS:
            nonzero = np.flatnonzero(values != 0)
            if nonzero.size == 0:
                continue
            starts = (nonzero + original_start).astype(np.int64)
            handle.addEntries(
                [contig] * len(starts),
                starts.tolist(),
                ends=(starts + 1).tolist(),
                values=values[nonzero].astype(float).tolist(),
            )
        else:
            handle.addEntries(
                contig,
                int(original_start),
                values=values.astype(float).tolist(),
                span=1,
                step=1,
            )


def _wig_value(track: str, value) -> str:
    if track in {
        "coverage",
        "fragment_ends",
        "fragment_left_ends",
        "fragment_right_ends",
    }:
        return str(int(value))
    return f"{float(value):.6g}"


def write_wig_tracks(
    scores,
    contig: str,
    adjusted_start: int,
    original_start: int,
    original_end: int,
    handles: Mapping,
    tracks: Sequence[str],
) -> None:
    for track in tracks:
        if track not in scores or track not in handles:
            continue
        handle = handles[track]
        array = scores[track][0][2]
        if track in SPARSE_TRACKS:
            handle.write(f"variableStep chrom={contig}\n")
            for position in range(original_start, original_end):
                index = position - adjusted_start
                if 0 <= index < len(array) and array[index] != 0:
                    handle.write(f"{position + 1}\t{_wig_value(track, array[index])}\n")
        else:
            handle.write(
                f"fixedStep chrom={contig} start={original_start + 1} step=1\n"
            )
            for position in range(original_start, original_end):
                index = position - adjusted_start
                value = array[index] if 0 <= index < len(array) else 0
                handle.write(_wig_value(track, value) + "\n")


def write_tracks(
    scores,
    contig: str,
    adjusted_start: int,
    original_start: int,
    original_end: int,
    tracks: Sequence[str],
    bigwig_handles,
    wig_handles,
) -> None:
    if bigwig_handles:
        write_bigwig_tracks(
            scores,
            contig,
            adjusted_start,
            original_start,
            original_end,
            bigwig_handles,
            tracks,
        )
    if wig_handles:
        write_wig_tracks(
            scores,
            contig,
            adjusted_start,
            original_start,
            original_end,
            wig_handles,
            tracks,
        )


def _write_completion_marker(final_path: Path) -> None:
    marker = final_path.with_name(final_path.name + ".complete.json")
    temporary = marker.with_name(marker.name + ".partial")
    temporary.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "complete": True,
                "output": str(final_path.resolve()),
                "output_size_bytes": int(final_path.stat().st_size),
                "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, marker)


def close_track_handles(*handle_groups, commit: bool = True) -> None:
    """Close writers and atomically publish every track only after success."""
    first_error: Exception | None = None
    for group in handle_groups:
        for handle in group.values():
            try:
                handle.close()
            except Exception as error:
                first_error = first_error or error
    if commit and first_error is None:
        for group in handle_groups:
            for partial_path, final_path in getattr(group, "paths", {}).values():
                os.replace(partial_path, final_path)
                _write_completion_marker(final_path)
    if first_error is not None and commit:
        raise first_error
