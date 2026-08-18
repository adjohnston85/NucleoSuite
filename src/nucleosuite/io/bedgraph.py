"""Validated streaming bedGraph writers used by suite staging mode."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

from nucleosuite.scoring.basic_tracks import SPARSE_TRACKS


STAGED_BEDGRAPH_SCHEMA_VERSION = 1


def staged_bedgraph_path(
    *,
    bigwig_path: str | Path,
    staging_root: str | Path,
    source_root: str | Path,
    source_id: str,
) -> Path:
    """Return the stable staged bedGraph path corresponding to one worker BigWig."""

    source = Path(bigwig_path).resolve()
    root = Path(source_root).resolve()
    try:
        relative = source.relative_to(root)
    except ValueError as error:
        raise ValueError(
            f"Track output {source} is not located beneath staged-bedGraph source root {root}"
        ) from error
    safe_id = str(source_id).replace("/", "_").replace("\\", "_").replace(":", "_")
    return Path(staging_root).resolve() / "per_contig" / safe_id / relative.with_suffix(".bedGraph")


@dataclass
class _PendingRun:
    chrom: str
    start: int
    end: int
    value: float


class ValidatedBedGraphWriter:
    """Write sorted, non-overlapping bedGraph records with run-length compression.

    Ordering is checked while records are generated, so no separate validation scan
    is required.  The final bedGraph is atomically published only after the writer
    closes successfully, followed by a small completion metadata file.
    """

    def __init__(
        self,
        path: str | Path,
        *,
        track: str,
        chrom_order: Sequence[str],
        source_bigwig: str | Path,
    ) -> None:
        self.path = Path(path)
        self.partial_path = self.path.with_name(self.path.name + ".partial")
        self.metadata_path = self.path.with_name(self.path.name + ".complete.json")
        self.track = str(track)
        self.source_bigwig = str(Path(source_bigwig).resolve())
        self.chrom_rank = {chrom: index for index, chrom in enumerate(chrom_order)}
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.partial_path.unlink(missing_ok=True)
        self.metadata_path.unlink(missing_ok=True)
        self.handle = self.partial_path.open("w", encoding="utf-8", newline="")
        self.pending: _PendingRun | None = None
        self.records = 0
        self.first_chrom: str | None = None
        self.first_start: int | None = None
        self.last_chrom: str | None = None
        self.last_start: int | None = None
        self.last_end: int | None = None
        self.closed = False

    def _validate_record(self, chrom: str, start: int, end: int) -> None:
        if chrom not in self.chrom_rank:
            raise ValueError(f"Unknown bedGraph chromosome {chrom!r} for {self.path}")
        if start < 0 or end <= start:
            raise ValueError(f"Invalid bedGraph interval {chrom}:{start}-{end} for {self.path}")
        if self.last_chrom is None:
            return
        previous_rank = self.chrom_rank[self.last_chrom]
        current_rank = self.chrom_rank[chrom]
        if current_rank < previous_rank:
            raise ValueError(
                f"bedGraph chromosome order inversion in {self.path}: "
                f"{chrom} follows {self.last_chrom}"
            )
        if chrom == self.last_chrom:
            if self.last_start is not None and start < self.last_start:
                raise ValueError(
                    f"bedGraph coordinate inversion in {self.path}: "
                    f"{chrom}:{start} follows {self.last_start}"
                )
            if self.last_end is not None and start < self.last_end:
                raise ValueError(
                    f"Overlapping bedGraph intervals in {self.path}: "
                    f"{chrom}:{start}-{end} starts before previous end {self.last_end}"
                )

    def _emit_pending(self) -> None:
        pending = self.pending
        if pending is None:
            return
        self._validate_record(pending.chrom, pending.start, pending.end)
        self.handle.write(
            f"{pending.chrom}\t{pending.start}\t{pending.end}\t{pending.value:.12g}\n"
        )
        if self.first_chrom is None:
            self.first_chrom = pending.chrom
            self.first_start = pending.start
        self.last_chrom = pending.chrom
        self.last_start = pending.start
        self.last_end = pending.end
        self.records += 1
        self.pending = None

    def add_interval(self, chrom: str, start: int, end: int, value: float) -> None:
        start_i = int(start)
        end_i = int(end)
        value_f = float(value)
        if not np.isfinite(value_f):
            raise ValueError(
                f"Non-finite bedGraph value for {chrom}:{start_i}-{end_i} in {self.path}"
            )
        pending = self.pending
        if (
            pending is not None
            and pending.chrom == chrom
            and pending.end == start_i
            and pending.value == value_f
        ):
            pending.end = end_i
            return
        self._emit_pending()
        self.pending = _PendingRun(chrom, start_i, end_i, value_f)

    def add_values(
        self,
        *,
        chrom: str,
        start: int,
        values: np.ndarray,
        sparse: bool,
    ) -> None:
        array = np.asarray(values, dtype=np.float64)
        if array.size == 0:
            return
        if sparse:
            indexes = np.flatnonzero(array != 0)
            for index in indexes:
                position = int(start) + int(index)
                self.add_interval(chrom, position, position + 1, float(array[index]))
            return

        run_start = 0
        current = float(array[0])
        for index in range(1, int(array.size)):
            value = float(array[index])
            if value != current:
                self.add_interval(
                    chrom,
                    int(start) + run_start,
                    int(start) + index,
                    current,
                )
                run_start = index
                current = value
        self.add_interval(
            chrom,
            int(start) + run_start,
            int(start) + int(array.size),
            current,
        )

    def close(self, *, commit: bool = True) -> None:
        if self.closed:
            return
        try:
            self._emit_pending()
            self.handle.flush()
            os.fsync(self.handle.fileno())
            self.handle.close()
            if not commit:
                self.closed = True
                return
            os.replace(self.partial_path, self.path)
            metadata = {
                "schema_version": STAGED_BEDGRAPH_SCHEMA_VERSION,
                "complete": True,
                "sorted": True,
                "nonoverlapping": True,
                "track": self.track,
                "source_bigwig": self.source_bigwig,
                "bedgraph": str(self.path.resolve()),
                "records": self.records,
                "first_chrom": self.first_chrom,
                "first_start": self.first_start,
                "last_chrom": self.last_chrom,
                "last_start": self.last_start,
                "last_end": self.last_end,
                "size_bytes": self.path.stat().st_size,
            }
            temporary_metadata = self.metadata_path.with_name(
                self.metadata_path.name + ".partial"
            )
            temporary_metadata.write_text(
                json.dumps(metadata, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            os.replace(temporary_metadata, self.metadata_path)
            self.closed = True
        except Exception:
            try:
                self.handle.close()
            except Exception:
                pass
            raise


def open_staged_bedgraph_handles(
    *,
    output_prefix: str,
    tracks: Iterable[str],
    staging_root: str | Path | None,
    source_root: str | Path | None,
    source_id: str | None,
    chrom_order: Sequence[str],
) -> dict[str, ValidatedBedGraphWriter]:
    if not staging_root:
        return {}
    if not source_root or not source_id:
        raise ValueError(
            "Staged bedGraph output requires source root and source identifier"
        )
    handles: dict[str, ValidatedBedGraphWriter] = {}
    for track in dict.fromkeys(tracks):
        bigwig_path = Path(f"{output_prefix}_{track}.bw")
        bedgraph_path = staged_bedgraph_path(
            bigwig_path=bigwig_path,
            staging_root=staging_root,
            source_root=source_root,
            source_id=source_id,
        )
        handles[track] = ValidatedBedGraphWriter(
            bedgraph_path,
            track=track,
            chrom_order=chrom_order,
            source_bigwig=bigwig_path,
        )
    return handles


def write_staged_bedgraph_tracks(
    *,
    scores: Mapping,
    contig: str,
    adjusted_start: int,
    original_start: int,
    original_end: int,
    handles: Mapping[str, ValidatedBedGraphWriter],
    tracks: Sequence[str],
) -> None:
    for track in tracks:
        if track not in scores or track not in handles:
            continue
        array = np.asarray(scores[track][0][2], dtype=np.float64)
        left = original_start - adjusted_start
        right = original_end - adjusted_start
        values = array[left:right]
        handles[track].add_values(
            chrom=contig,
            start=original_start,
            values=values,
            sparse=track in SPARSE_TRACKS,
        )


def close_staged_bedgraph_handles(
    *groups: Mapping[str, ValidatedBedGraphWriter], commit: bool = True
) -> None:
    first_error: Exception | None = None
    for group in groups:
        for handle in group.values():
            try:
                handle.close(commit=commit)
            except Exception as error:  # preserve cleanup attempts for all tracks
                if first_error is None:
                    first_error = error
    if first_error is not None:
        raise first_error
