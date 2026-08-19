#!/usr/bin/env python3
"""Compare summit positions and scores between two BED interval sets."""

from __future__ import annotations

import argparse
import bisect
import csv
import gc
import heapq
import json
import math
import os
import sys
from array import array
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy import stats

from nucleosuite.io import open_text as open_interval_text
from nucleosuite.core.regions import canonical_contig_key
from nucleosuite.core.blacklist import BlacklistIndex, load_blacklist_unbounded
from nucleosuite.progress import ProgressReporter
from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter


@dataclass(frozen=True)
class PositionRecord:
    """One BED interval with a resolved absolute summit and numeric score."""

    source: str
    chrom: str
    start: int
    end: int
    summit: int
    score: float
    name: str
    line_number: int


@dataclass(frozen=True)
class MatchedPair:
    """One matched A/B position pair."""

    a: PositionRecord
    b: PositionRecord
    query_source: str

    @property
    def signed_distance(self) -> int:
        return self.b.summit - self.a.summit

    @property
    def absolute_distance(self) -> int:
        return abs(self.signed_distance)


@dataclass(frozen=True)
class MatchResult:
    """Matched records and matching diagnostics."""

    pairs: list[MatchedPair]
    query_source: str
    target_source: str
    query_count: int
    target_count: int
    unmatched_no_target_chrom: int
    unmatched_distance: int
    unmatched_unique: int


@dataclass(frozen=True)
class PercentileAssignment:
    """Independent score-percentile assignment for one input position."""

    percentile: int
    group_lower: int
    group_upper: int

    @property
    def label(self) -> str:
        return f"{self.group_lower}-{self.group_upper}"


@dataclass
class _CandidateCursor:
    """Generate target indices in increasing distance from one query summit."""

    query: PositionRecord
    targets: list[PositionRecord]
    positions: list[int]
    left: int
    right: int

    @classmethod
    def create(
        cls,
        query: PositionRecord,
        targets: list[PositionRecord],
        positions: list[int],
    ) -> "_CandidateCursor":
        right = bisect.bisect_left(positions, query.summit)
        return cls(query, targets, positions, right - 1, right)

    def next(self) -> tuple[int, int] | None:
        """Return ``(absolute_distance, target_index)`` for the next candidate."""

        if self.left < 0 and self.right >= len(self.targets):
            return None
        if self.left < 0:
            index = self.right
            self.right += 1
            return abs(self.positions[index] - self.query.summit), index
        if self.right >= len(self.targets):
            index = self.left
            self.left -= 1
            return abs(self.positions[index] - self.query.summit), index

        left_distance = abs(self.positions[self.left] - self.query.summit)
        right_distance = abs(self.positions[self.right] - self.query.summit)
        if left_distance <= right_distance:
            index = self.left
            self.left -= 1
            return left_distance, index
        index = self.right
        self.right += 1
        return right_distance, index


@dataclass
class _PositionNode:
    """One coordinate group in the ordered one-to-one matching frontier."""

    kind: str
    coordinate: int
    record_indices: deque[int]
    previous: int | None = None
    following: int | None = None
    version: int = 0
    alive: bool = True


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite compare-positions",
        description=(
            "Match the smaller of two BED position sets to nearest positions in "
            "the other set and quantify positional and score agreement."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("--bed-a", required=True, help="Method A BED, BED.gz, or bigBed file.")
    parser.add_argument("--bed-b", required=True, help="Method B BED, BED.gz, or bigBed file.")
    parser.add_argument(
        "--blacklist-bed",
        help="BED blacklist; complete overlapping records are excluded from both inputs.",
    )
    parser.add_argument(
        "--summit-column-a",
        type=int,
        default=None,
        help="One-based absolute summit column for BED A; default uses the BED midpoint.",
    )
    parser.add_argument(
        "--summit-column-b",
        type=int,
        default=None,
        help="One-based absolute summit column for BED B; default uses the BED midpoint.",
    )
    parser.add_argument(
        "--score-column-a",
        type=int,
        default=5,
        help="One-based numeric score column for BED A.",
    )
    parser.add_argument(
        "--score-column-b",
        type=int,
        default=5,
        help="One-based numeric score column for BED B.",
    )
    parser.add_argument("--label-a", default=None, help="Display label for method A.")
    parser.add_argument("--label-b", default=None, help="Display label for method B.")
    parser.add_argument(
        "--matching",
        choices=("one-to-one", "many-to-one", "unique"),
        default="one-to-one",
        help=(
            "Nearest-neighbour matching policy. 'one-to-one' uses distance-prioritized "
            "greedy matching without target reuse; 'many-to-one' permits reuse. "
            "'unique' is an accepted alias for one-to-one."
        ),
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=None,
        help="Maximum allowed absolute summit distance in bp; omit for no limit.",
    )
    parser.add_argument(
        "--distance-bins",
        default="5,10,20,50,100",
        help=(
            "Comma-separated upper bounds for distance-bin summaries. A final "
            "open-ended bin is added automatically."
        ),
    )
    parser.add_argument(
        "--plot-max-points",
        type=int,
        default=200000,
        help="Maximum matched pairs drawn in scatter plots; 0 draws all pairs.",
    )
    parser.add_argument("--plot-seed", type=int, default=1, help="Scatter-plot subsampling seed.")
    parser.add_argument(
        "--score-normalization",
        choices=("raw", "zscore", "percentile"),
        default="zscore",
        help="Score representation used for plots, correlations, regression, and distance-bin statistics.",
    )
    parser.add_argument(
        "--correlation-method",
        choices=("spearman", "pearson", "both"),
        default="spearman",
        help="Correlation statistic shown in the distance-bin plot and emphasized in summaries.",
    )
    parser.add_argument(
        "--score-z-limit",
        type=float,
        default=10.0,
        help=(
            "Symmetric absolute x- and y-axis limit for the score-correlation plot "
            "when --score-normalization zscore is used. Use 0 to disable the limit."
        ),
    )
    parser.add_argument(
        "--histogram-bin-width",
        type=float,
        default=1.0,
        help="Distance-histogram bin width in bp.",
    )
    parser.add_argument(
        "--histogram-x-min",
        type=float,
        default=0.0,
        help="Displayed lower x-axis limit for the distance histogram in bp.",
    )
    parser.add_argument(
        "--histogram-x-max",
        type=float,
        default=300.0,
        help="Displayed upper x-axis limit for the distance histogram in bp.",
    )
    parser.add_argument(
        "--distance-x-major-tick",
        type=float,
        default=None,
        help="Major x-axis tick interval in bp for numeric distance plots; default is automatic.",
    )
    parser.add_argument(
        "--distance-x-minor-tick",
        type=float,
        default=None,
        help=(
            "Minor x-axis tick interval in bp. Automatic defaults use 5 bp for a "
            "10 bp major interval, 10 bp when the major interval is greater than 50 bp, "
            "and half the major interval otherwise."
        ),
    )
    parser.add_argument(
        "--percentile-interval",
        type=int,
        default=25,
        help=(
            "Width of independently calculated score-percentile groups. The default "
            "quartiles are 0-25, 25-50, 50-75, and 75-100. Groups use lower-exclusive, "
            "upper-inclusive boundaries except at zero. Two directional analyses are "
            "produced: each A percentile group versus all B positions, and each B "
            "percentile group versus all A positions."
        ),
    )
    parser.add_argument(
        "--percentile-boxplot-y-max",
        type=float,
        default=500.0,
        help=(
            "Displayed upper y-axis limit in bp for directional percentile-distance "
            "boxplots. Use 0 to disable the limit."
        ),
    )
    parser.add_argument(
        "--skip-percentile-distance-analysis",
        action="store_true",
        help=(
            "Skip the two directional score-percentile distance analyses and their "
            "box-and-whisker plots."
        ),
    )
    parser.add_argument(
        "--skip-pairs-tsv",
        action="store_true",
        help="Run the main analyses without writing the potentially large main pair TSV.",
    )
    parser.add_argument(
        "--skip-percentile-pairs-tsv",
        action="store_true",
        help=(
            "Run directional percentile summaries and plots without writing their "
            "potentially large detailed pair TSVs."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Ignore compatible stage checkpoints and recompute every analysis stage.",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Figure resolution setting; --plot-dpi takes precedence when supplied.")
    parser.add_argument(
        "-o",
        "--output-prefix",
        default=None,
        help="Output prefix. Default combines the two input basenames.",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress progress messages.")
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(
        parser,
        cores_option="--memory-intensive-analysis-cores",
        cores_help=(
            "Concurrent contig workers for this memory-intensive analysis. This "
            "budget is independent of suite --cores (default: 1)."
        ),
    )
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def _positive_column(value: int | None, option: str) -> int | None:
    if value is not None and value < 1:
        raise ValueError(f"{option} must be a one-based column number of at least 1.")
    return value


def _parse_integer_coordinate(value: str, context: str) -> int:
    try:
        numeric = float(value)
    except ValueError as exc:
        raise ValueError(f"{context}: expected a numeric summit coordinate, found {value!r}.") from exc
    if not math.isfinite(numeric) or not numeric.is_integer():
        raise ValueError(f"{context}: summit coordinates must be finite integers, found {value!r}.")
    return int(numeric)


def read_positions(
    path: str | Path,
    source: str,
    summit_column: int | None,
    score_column: int,
    blacklist: BlacklistIndex | None = None,
    excluded_counter: list[int] | None = None,
    progress: ProgressReporter | None = None,
) -> list[PositionRecord]:
    """Read positions from BED, BED.gz, or bigBed using one-based column options."""

    input_path = Path(path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input interval file not found: {input_path}")
    _positive_column(summit_column, f"--summit-column-{source.lower()}")
    _positive_column(score_column, f"--score-column-{source.lower()}")

    required_column = max(3, score_column, summit_column or 0)
    records: list[PositionRecord] = []
    if progress is not None:
        progress.file_start(source, input_path)
    seen_contigs: set[str] = set()
    with open_interval_text(input_path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            text = raw.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t") if "\t" in text else text.split()
            if len(fields) < required_column:
                raise ValueError(
                    f"{input_path}:{line_number}: requires at least {required_column} "
                    f"columns for the selected summit and score fields; found {len(fields)}."
                )
            chrom = fields[0]
            if progress is not None and chrom not in seen_contigs:
                seen_contigs.add(chrom)
                progress.reading_contig(source, chrom)
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{input_path}:{line_number}: BED start and end must be integers."
                ) from exc
            if start < 0 or end <= start:
                raise ValueError(
                    f"{input_path}:{line_number}: require 0 <= start < end."
                )
            if blacklist is not None and blacklist.overlaps(chrom, start, end):
                if excluded_counter is not None:
                    excluded_counter[0] += 1
                continue
            if summit_column is None:
                summit = (start + end) // 2
            else:
                summit = _parse_integer_coordinate(
                    fields[summit_column - 1],
                    f"{input_path}:{line_number}",
                )
            try:
                score = float(fields[score_column - 1])
            except ValueError as exc:
                raise ValueError(
                    f"{input_path}:{line_number}: score column {score_column} "
                    f"must be numeric."
                ) from exc
            if not math.isfinite(score):
                raise ValueError(
                    f"{input_path}:{line_number}: score column {score_column} must be finite."
                )
            name = fields[3] if len(fields) >= 4 else f"{chrom}:{start}-{end}"
            records.append(
                PositionRecord(
                    source=source,
                    chrom=chrom,
                    start=start,
                    end=end,
                    summit=summit,
                    score=score,
                    name=name,
                    line_number=line_number,
                )
            )

    if not records:
        raise ValueError(f"No interval records were read from {input_path}.")
    if progress is not None:
        progress.file_complete(source, input_path, len(records))
    return records


def _records_by_chrom(records: Iterable[PositionRecord]) -> dict[str, list[PositionRecord]]:
    grouped: dict[str, list[PositionRecord]] = defaultdict(list)
    for record in records:
        grouped[canonical_contig_key(record.chrom)].append(record)
    for chrom in grouped:
        grouped[chrom].sort(
            key=lambda record: (
                record.summit,
                record.start,
                record.end,
                record.line_number,
            )
        )
    return dict(grouped)


def _nearest_index(query: PositionRecord, targets: list[PositionRecord], positions: list[int]) -> int:
    insertion = bisect.bisect_left(positions, query.summit)
    candidates: list[int] = []
    if insertion > 0:
        candidates.append(insertion - 1)
    if insertion < len(targets):
        candidates.append(insertion)
    return min(
        candidates,
        key=lambda index: (
            abs(targets[index].summit - query.summit),
            targets[index].summit,
            targets[index].start,
            targets[index].end,
            targets[index].line_number,
        ),
    )


def _make_pair(query: PositionRecord, target: PositionRecord, query_source: str) -> MatchedPair:
    if query_source == "A":
        return MatchedPair(a=query, b=target, query_source="A")
    return MatchedPair(a=target, b=query, query_source="B")


def match_many_to_one(
    query_records: list[PositionRecord],
    target_records: list[PositionRecord],
    query_source: str,
    max_distance: float | None,
    progress: ProgressReporter | None = None,
    progress_stage: str = "Matching many-to-one",
) -> MatchResult:
    target_by_chrom = _records_by_chrom(target_records)
    target_positions = {
        chrom: [record.summit for record in records]
        for chrom, records in target_by_chrom.items()
    }
    pairs: list[MatchedPair] = []
    unmatched_no_target_chrom = 0
    unmatched_distance = 0

    ordered_queries = sorted(
        query_records,
        key=lambda record: (
            canonical_contig_key(record.chrom),
            record.summit,
            record.line_number,
        ),
    )
    last_chrom = None
    contig_index = 0
    total_contigs = len({canonical_contig_key(record.chrom) for record in ordered_queries})
    for query in ordered_queries:
        query_chrom_key = canonical_contig_key(query.chrom)
        if progress is not None and query_chrom_key != last_chrom:
            contig_index += 1
            last_chrom = query_chrom_key
            progress.contig(
                progress_stage, query.chrom, contig_index, total_contigs
            )
        targets = target_by_chrom.get(query_chrom_key)
        if not targets:
            unmatched_no_target_chrom += 1
            continue
        target_index = _nearest_index(query, targets, target_positions[query_chrom_key])
        target = targets[target_index]
        distance = abs(target.summit - query.summit)
        if max_distance is not None and distance > max_distance:
            unmatched_distance += 1
            continue
        pairs.append(_make_pair(query, target, query_source))

    return MatchResult(
        pairs=pairs,
        query_source=query_source,
        target_source="B" if query_source == "A" else "A",
        query_count=len(query_records),
        target_count=len(target_records),
        unmatched_no_target_chrom=unmatched_no_target_chrom,
        unmatched_distance=unmatched_distance,
        unmatched_unique=0,
    )


def match_unique(
    query_records: list[PositionRecord],
    target_records: list[PositionRecord],
    query_source: str,
    max_distance: float | None,
    progress: ProgressReporter | None = None,
    progress_stage: str = "Matching one-to-one",
) -> MatchResult:
    """Distance-prioritized greedy one-to-one matching in near-linear space.

    After exact-coordinate pairs are resolved, the closest remaining
    query/target pair must lie on a boundary between adjacent opposite-type
    coordinate groups. A heap stores only those boundaries. Removing a matched
    pair exposes at most two new boundaries, giving approximately
    ``O((Q + T) log(Q + T))`` time instead of repeatedly scanning targets that
    have already been assigned.
    """

    query_by_chrom = _records_by_chrom(query_records)
    target_by_chrom = _records_by_chrom(target_records)
    pairs: list[MatchedPair] = []
    unmatched_no_target_chrom = 0
    unmatched_distance = 0
    unmatched_unique = 0

    contigs = list(query_by_chrom)
    for chrom_index, (chrom, queries) in enumerate(query_by_chrom.items(), start=1):
        if progress is not None:
            progress.contig(
                progress_stage,
                queries[0].chrom if queries else chrom,
                chrom_index,
                len(contigs),
                len(queries),
            )
        targets = target_by_chrom.get(chrom)
        if not targets:
            unmatched_no_target_chrom += len(queries)
            continue

        queries_at: dict[int, deque[int]] = defaultdict(deque)
        targets_at: dict[int, deque[int]] = defaultdict(deque)
        for query_index, record in enumerate(queries):
            queries_at[record.summit].append(query_index)
        for target_index, record in enumerate(targets):
            targets_at[record.summit].append(target_index)

        matched_query_count = 0
        nodes: list[_PositionNode] = []
        for coordinate in sorted(set(queries_at) | set(targets_at)):
            query_indexes = queries_at.get(coordinate, deque())
            target_indexes = targets_at.get(coordinate, deque())
            while query_indexes and target_indexes:
                query_index = query_indexes.popleft()
                target_index = target_indexes.popleft()
                pairs.append(
                    _make_pair(
                        queries[query_index], targets[target_index], query_source
                    )
                )
                matched_query_count += 1
            if query_indexes:
                nodes.append(_PositionNode("Q", coordinate, query_indexes))
            elif target_indexes:
                nodes.append(_PositionNode("T", coordinate, target_indexes))

        for node_index, node in enumerate(nodes):
            node.previous = node_index - 1 if node_index else None
            node.following = node_index + 1 if node_index + 1 < len(nodes) else None

        boundary_heap: list[tuple[int, int, int, int, int, int, int]] = []

        def push_boundary(left_index: int | None) -> None:
            if left_index is None:
                return
            left = nodes[left_index]
            right_index = left.following
            if not left.alive or right_index is None:
                return
            right = nodes[right_index]
            if not right.alive or left.kind == right.kind:
                return
            if left.kind == "Q":
                query_index = left.record_indices[0]
                target_index = right.record_indices[0]
            else:
                query_index = right.record_indices[0]
                # The bisect cursor approaches a duplicate target
                # coordinate from the right by walking leftward, so it selects
                # the last target at that coordinate first.
                target_index = left.record_indices[-1]
            heapq.heappush(
                boundary_heap,
                (
                    abs(right.coordinate - left.coordinate),
                    query_index,
                    target_index,
                    left_index,
                    right_index,
                    left.version,
                    right.version,
                ),
            )

        for node_index in range(max(0, len(nodes) - 1)):
            push_boundary(node_index)

        def unlink(node_index: int) -> None:
            node = nodes[node_index]
            previous = node.previous
            following = node.following
            if previous is not None:
                nodes[previous].following = following
            if following is not None:
                nodes[following].previous = previous
            node.alive = False
            node.previous = None
            node.following = None
            node.version += 1

        while boundary_heap:
            (
                distance,
                query_index,
                target_index,
                left_index,
                right_index,
                left_version,
                right_version,
            ) = heapq.heappop(boundary_heap)
            left = nodes[left_index]
            right = nodes[right_index]
            if (
                not left.alive
                or not right.alive
                or left.following != right_index
                or right.previous != left_index
                or left.version != left_version
                or right.version != right_version
            ):
                continue
            if max_distance is not None and distance > max_distance:
                break
            current_query = (
                left.record_indices[0]
                if left.kind == "Q"
                else right.record_indices[0]
            )
            current_target = (
                right.record_indices[0]
                if left.kind == "Q"
                else left.record_indices[-1]
            )
            if current_query != query_index or current_target != target_index:
                continue

            pairs.append(
                _make_pair(
                    queries[query_index], targets[target_index], query_source
                )
            )
            matched_query_count += 1
            affected = {
                left.previous,
                left_index,
                right_index,
                right.following,
            }
            if left.kind == "Q":
                left.record_indices.popleft()
                right.record_indices.popleft()
            else:
                left.record_indices.pop()
                right.record_indices.popleft()
            left.version += 1
            right.version += 1
            if not left.record_indices:
                unlink(left_index)
            if not right.record_indices:
                unlink(right_index)
            for candidate_index in list(affected):
                if candidate_index is None or not nodes[candidate_index].alive:
                    continue
                push_boundary(nodes[candidate_index].previous)
                push_boundary(candidate_index)

        remaining_query_indices: list[int] = []
        for node in nodes:
            if not node.alive:
                continue
            if node.kind == "Q":
                remaining_query_indices.extend(node.record_indices)
        if max_distance is None:
            unmatched_unique += len(remaining_query_indices)
        else:
            # Preserve the former candidate-cursor diagnostic semantics: an
            # unmatched query kept expanding through consumed targets.  It was
            # classified as distance-rejected when that expansion first reached
            # any candidate beyond the limit, and as a pure uniqueness conflict
            # only when every original target was within the limit.
            first_target = targets[0].summit
            last_target = targets[-1].summit
            for query_index in remaining_query_indices:
                summit = queries[query_index].summit
                farthest = max(
                    abs(first_target - summit), abs(last_target - summit)
                )
                if farthest > max_distance:
                    unmatched_distance += 1
                else:
                    unmatched_unique += 1

        expected_unmatched = len(queries) - matched_query_count
        observed_unmatched = (
            len(remaining_query_indices)
        )
        if expected_unmatched != observed_unmatched:
            raise RuntimeError(
                "Internal one-to-one matching count mismatch: "
                f"expected {expected_unmatched}, observed {observed_unmatched}"
            )

    return MatchResult(
        pairs=sorted(
            pairs,
            key=lambda pair: (
                pair.a.chrom,
                pair.a.summit if query_source == "A" else pair.b.summit,
                pair.a.line_number,
                pair.b.line_number,
            ),
        ),
        query_source=query_source,
        target_source="B" if query_source == "A" else "A",
        query_count=len(query_records),
        target_count=len(target_records),
        unmatched_no_target_chrom=unmatched_no_target_chrom,
        unmatched_distance=unmatched_distance,
        unmatched_unique=unmatched_unique,
    )


def match_positions(
    records_a: list[PositionRecord],
    records_b: list[PositionRecord],
    matching: str,
    max_distance: float | None,
    progress: ProgressReporter | None = None,
    progress_stage: str | None = None,
) -> MatchResult:
    """Use the smaller valid record set as query and match by chromosome."""

    if max_distance is not None and max_distance < 0:
        raise ValueError("--max-distance must be zero or greater.")
    if len(records_a) <= len(records_b):
        query_records, target_records, query_source = records_a, records_b, "A"
    else:
        query_records, target_records, query_source = records_b, records_a, "B"

    if matching == "many-to-one":
        return match_many_to_one(
            query_records,
            target_records,
            query_source,
            max_distance,
            progress=progress,
            progress_stage=progress_stage or "Matching many-to-one",
        )
    if matching in {"unique", "one-to-one"}:
        return match_unique(
            query_records,
            target_records,
            query_source,
            max_distance,
            progress=progress,
            progress_stage=progress_stage or "Matching one-to-one",
        )
    raise ValueError(f"Unknown matching policy: {matching}")


def match_positions_directional(
    query_records: Sequence[PositionRecord],
    target_records: Sequence[PositionRecord],
    query_source: str,
    matching: str,
    max_distance: float | None,
    progress: ProgressReporter | None = None,
    progress_stage: str | None = None,
) -> MatchResult:
    """Match a designated query set against a designated complete target set.

    Unlike :func:`match_positions`, this function never swaps the two inputs based
    on their sizes. It is used by the directional percentile analyses, where the
    percentile group must remain the query and the complete opposite callset must
    remain the target.
    """

    if max_distance is not None and max_distance < 0:
        raise ValueError("--max-distance must be zero or greater.")
    if query_source not in {"A", "B"}:
        raise ValueError("query_source must be 'A' or 'B'.")

    query_list = list(query_records)
    target_list = list(target_records)
    if matching == "many-to-one":
        return match_many_to_one(
            query_list,
            target_list,
            query_source,
            max_distance,
            progress=progress,
            progress_stage=progress_stage or "Matching many-to-one",
        )
    if matching in {"unique", "one-to-one"}:
        return match_unique(
            query_list,
            target_list,
            query_source,
            max_distance,
            progress=progress,
            progress_stage=progress_stage or "Matching one-to-one",
        )
    raise ValueError(f"Unknown matching policy: {matching}")


def _validate_percentile_interval(interval: int) -> int:
    if interval < 1 or interval > 100:
        raise ValueError("--percentile-interval must be between 1 and 100.")
    return interval


def _percentile_group_bounds(interval: int) -> list[tuple[int, int, str]]:
    """Return boundary-labelled groups such as 0-25, 25-50, 50-75, 75-100.

    Membership is lower-exclusive and upper-inclusive, except that the first
    boundary begins at zero. Score percentiles themselves remain integer values
    from 1 through 100.
    """

    _validate_percentile_interval(interval)
    groups: list[tuple[int, int, str]] = []
    lower = 0
    while lower < 100:
        upper = min(100, lower + interval)
        groups.append((lower, upper, f"{lower}-{upper}"))
        lower = upper
    return groups


def assign_score_percentiles(
    records: Sequence[PositionRecord],
    interval: int,
) -> dict[int, PercentileAssignment]:
    """Assign equal-frequency score percentiles independently within one BED.

    Scores are ordered from low to high. Ties are resolved deterministically by
    genomic position and input line number so groups remain close to equal size.
    The returned mapping is keyed by the input line number, which is unique within
    each BED file.
    """

    _validate_percentile_interval(interval)
    ordered = sorted(
        records,
        key=lambda record: (
            record.score,
            record.chrom,
            record.summit,
            record.start,
            record.end,
            record.line_number,
        ),
    )
    count = len(ordered)
    assignments: dict[int, PercentileAssignment] = {}
    for rank, record in enumerate(ordered, start=1):
        percentile = int(math.ceil(rank * 100.0 / count))
        lower = ((percentile - 1) // interval) * interval
        upper = min(100, lower + interval)
        assignments[record.line_number] = PercentileAssignment(
            percentile=percentile,
            group_lower=lower,
            group_upper=upper,
        )
    return assignments


def assign_score_percentiles_compact(
    records: Sequence[PositionRecord],
    interval: int,
) -> np.ndarray:
    """Assign 1-100 score percentiles in a compact line-number-indexed array.

    Numeric lexicographic sorting preserves the public tie rule (genomic
    position followed by input line number) without constructing millions of
    Python mapping entries.
    """

    _validate_percentile_interval(interval)
    if not records:
        return np.zeros(0, dtype=np.uint8)
    count = len(records)
    scores = np.fromiter(
        (record.score for record in records), dtype=np.float64, count=count
    )
    chrom_codes = {chrom: index for index, chrom in enumerate(
        sorted({record.chrom for record in records})
    )}
    chroms = np.fromiter(
        (chrom_codes[record.chrom] for record in records), dtype=np.int32, count=count
    )
    summits = np.fromiter(
        (record.summit for record in records), dtype=np.int64, count=count
    )
    starts = np.fromiter(
        (record.start for record in records), dtype=np.int64, count=count
    )
    ends = np.fromiter(
        (record.end for record in records), dtype=np.int64, count=count
    )
    line_numbers = np.fromiter(
        (record.line_number for record in records), dtype=np.int64, count=count
    )
    # Use deterministic ordering by score, chromosome,
    # summit, interval start/end, then input line number.
    order = np.lexsort((line_numbers, ends, starts, summits, chroms, scores))
    percentiles_by_record = np.empty(len(records), dtype=np.uint8)
    ranks = np.arange(1, len(records) + 1, dtype=np.float64)
    percentiles_by_record[order] = np.ceil(
        ranks * 100.0 / len(records)
    ).astype(np.uint8)
    assignments = np.zeros(
        max(record.line_number for record in records) + 1, dtype=np.uint8
    )
    assignments[line_numbers] = percentiles_by_record
    return assignments


def directional_percentile_distance_rows(
    records_a: Sequence[PositionRecord],
    records_b: Sequence[PositionRecord],
    percentile_source: str,
    matching: str,
    max_distance: float | None,
    interval: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, np.ndarray]]:
    """Compare each percentile group from one source with all opposite positions.

    ``percentile_source='A'`` produces A percentile groups versus the complete B
    callset. ``percentile_source='B'`` produces the reciprocal analysis. The
    opposite callset is never split or filtered by percentile for matching.
    """

    if percentile_source not in {"A", "B"}:
        raise ValueError("percentile_source must be 'A' or 'B'.")

    groups = _percentile_group_bounds(interval)
    assignment_a = assign_score_percentiles(records_a, interval)
    assignment_b = assign_score_percentiles(records_b, interval)
    source_records = list(records_a if percentile_source == "A" else records_b)
    target_records = list(records_b if percentile_source == "A" else records_a)
    source_assignments = assignment_a if percentile_source == "A" else assignment_b
    target_source = "B" if percentile_source == "A" else "A"

    grouped_source: dict[str, list[PositionRecord]] = defaultdict(list)
    for record in source_records:
        grouped_source[source_assignments[record.line_number].label].append(record)

    detail_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    distances_by_group: dict[str, np.ndarray] = {}
    pair_number = 0

    for lower, upper, label in groups:
        query_group = grouped_source.get(label, [])
        pairs: list[MatchedPair] = []
        unmatched_no_target = 0
        unmatched_distance = 0
        unmatched_unique = 0
        if query_group:
            result = match_positions_directional(
                query_group,
                target_records,
                percentile_source,
                matching,
                max_distance,
            )
            pairs = result.pairs
            unmatched_no_target = result.unmatched_no_target_chrom
            unmatched_distance = result.unmatched_distance
            unmatched_unique = result.unmatched_unique

        matched_source = {
            pair.a.line_number if percentile_source == "A" else pair.b.line_number
            for pair in pairs
        }
        matched_target = {
            pair.b.line_number if percentile_source == "A" else pair.a.line_number
            for pair in pairs
        }
        distances = np.asarray([pair.absolute_distance for pair in pairs], dtype=np.float64)
        distances_by_group[label] = distances

        for pair in pairs:
            pair_number += 1
            a_assignment = assignment_a[pair.a.line_number]
            b_assignment = assignment_b[pair.b.line_number]
            signed_target_minus_query = (
                pair.signed_distance
                if percentile_source == "A"
                else -pair.signed_distance
            )
            detail_rows.append(
                {
                    "pair_id": f"percentile_pair_{pair_number:07d}",
                    "analysis_direction": f"{percentile_source}_percentiles_vs_all_{target_source}",
                    "percentile_source": percentile_source,
                    "target_source": target_source,
                    "percentile_group": label,
                    "group_lower_percentile": lower,
                    "group_upper_percentile": upper,
                    "query_source": percentile_source,
                    "chrom": pair.a.chrom,
                    "a_start": pair.a.start,
                    "a_end": pair.a.end,
                    "a_name": pair.a.name,
                    "a_summit": pair.a.summit,
                    "a_score": pair.a.score,
                    "a_score_percentile": a_assignment.percentile,
                    "b_start": pair.b.start,
                    "b_end": pair.b.end,
                    "b_name": pair.b.name,
                    "b_summit": pair.b.summit,
                    "b_score": pair.b.score,
                    "b_score_percentile": b_assignment.percentile,
                    "signed_distance_b_minus_a": pair.signed_distance,
                    "signed_distance_target_minus_query": signed_target_minus_query,
                    "absolute_distance": pair.absolute_distance,
                }
            )

        if distances.size:
            q1, median, q3 = np.percentile(distances, [25, 50, 75])
            minimum = float(np.min(distances))
            maximum = float(np.max(distances))
            mean = float(np.mean(distances))
            standard_deviation = float(np.std(distances, ddof=0))
        else:
            minimum = q1 = median = mean = q3 = maximum = standard_deviation = float("nan")

        summary_rows.append(
            {
                "analysis_direction": f"{percentile_source}_percentiles_vs_all_{target_source}",
                "percentile_source": percentile_source,
                "target_source": target_source,
                "percentile_group": label,
                "group_lower_percentile": lower,
                "group_upper_percentile": upper,
                "percentile_source_position_count": len(query_group),
                "all_target_position_count": len(target_records),
                "matched_pair_count": len(pairs),
                "matched_unique_source_positions": len(matched_source),
                "matched_unique_target_positions": len(matched_target),
                "unmatched_source_positions": len(query_group) - len(matched_source),
                "target_positions_not_used": len(target_records) - len(matched_target),
                "unmatched_query_no_target_chromosome": unmatched_no_target,
                "unmatched_query_beyond_maximum_distance": unmatched_distance,
                "unmatched_query_unique_assignment": unmatched_unique,
                "minimum_absolute_distance": minimum,
                "q1_absolute_distance": q1,
                "median_absolute_distance": median,
                "mean_absolute_distance": mean,
                "q3_absolute_distance": q3,
                "maximum_absolute_distance": maximum,
                "absolute_distance_standard_deviation": standard_deviation,
            }
        )

    return detail_rows, summary_rows, distances_by_group


PERCENTILE_DETAIL_FIELDS = [
    "pair_id",
    "analysis_direction",
    "percentile_source",
    "target_source",
    "percentile_group",
    "group_lower_percentile",
    "group_upper_percentile",
    "query_source",
    "chrom",
    "a_start",
    "a_end",
    "a_name",
    "a_summit",
    "a_score",
    "a_score_percentile",
    "b_start",
    "b_end",
    "b_name",
    "b_summit",
    "b_score",
    "b_score_percentile",
    "signed_distance_b_minus_a",
    "signed_distance_target_minus_query",
    "absolute_distance",
]

PERCENTILE_SUMMARY_FIELDS = [
    "analysis_direction",
    "percentile_source",
    "target_source",
    "percentile_group",
    "group_lower_percentile",
    "group_upper_percentile",
    "percentile_source_position_count",
    "all_target_position_count",
    "matched_pair_count",
    "matched_unique_source_positions",
    "matched_unique_target_positions",
    "unmatched_source_positions",
    "target_positions_not_used",
    "unmatched_query_no_target_chromosome",
    "unmatched_query_beyond_maximum_distance",
    "unmatched_query_unique_assignment",
    "minimum_absolute_distance",
    "q1_absolute_distance",
    "median_absolute_distance",
    "mean_absolute_distance",
    "q3_absolute_distance",
    "maximum_absolute_distance",
    "absolute_distance_standard_deviation",
]


def _atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def _save_numpy_atomic(path: Path, values: np.ndarray) -> None:
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("wb") as handle:
        np.save(handle, values, allow_pickle=False)
    os.replace(temporary, path)


def _file_signature(path: str | Path) -> dict[str, object]:
    value = Path(path).resolve()
    stat = value.stat()
    return {
        "path": str(value),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def _comparison_signature(args: argparse.Namespace) -> str:
    payload = {
        "schema_version": 1,
        "bed_a": _file_signature(args.bed_a),
        "bed_b": _file_signature(args.bed_b),
        "blacklist": (
            _file_signature(args.blacklist_bed)
            if getattr(args, "blacklist_bed", None)
            else None
        ),
        "summit_column_a": args.summit_column_a,
        "summit_column_b": args.summit_column_b,
        "score_column_a": args.score_column_a,
        "score_column_b": args.score_column_b,
        "label_a": args.label_a,
        "label_b": args.label_b,
        "matching": "one-to-one" if args.matching == "unique" else args.matching,
        "max_distance": args.max_distance,
        "percentile_interval": int(getattr(args, "percentile_interval", 25)),
        "distance_bins": args.distance_bins,
        "score_normalization": args.score_normalization,
        "correlation_method": args.correlation_method,
        "histogram_bin_width": args.histogram_bin_width,
        "histogram_x_min": args.histogram_x_min,
        "histogram_x_max": args.histogram_x_max,
        "score_z_limit": getattr(args, "score_z_limit", 10.0),
        "distance_x_major_tick": getattr(args, "distance_x_major_tick", None),
        "distance_x_minor_tick": getattr(args, "distance_x_minor_tick", None),
        "plot_max_points": args.plot_max_points,
        "plot_seed": args.plot_seed,
        "dpi": args.dpi,
        "percentile_boxplot_y_max": getattr(
            args, "percentile_boxplot_y_max", 500.0
        ),
        "skip_pairs_tsv": bool(getattr(args, "skip_pairs_tsv", False)),
        "skip_percentile_pairs_tsv": bool(
            getattr(args, "skip_percentile_pairs_tsv", False)
        ),
        "skip_percentile_distance_analysis": bool(
            getattr(args, "skip_percentile_distance_analysis", False)
        ),
    }
    return __import__("hashlib").sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _checkpoint_matches(
    marker: Path,
    *,
    signature: str,
    required_paths: Sequence[Path],
) -> bool:
    if not marker.is_file():
        return False
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
        if not payload.get("complete") or payload.get("signature") != signature:
            return False
        sizes = payload.get("sizes", {})
        for path in required_paths:
            if not path.is_file() or int(sizes.get(path.name, -1)) != path.stat().st_size:
                return False
        return True
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return False


def _write_checkpoint(
    marker: Path,
    *,
    signature: str,
    required_paths: Sequence[Path],
    stage: str,
) -> None:
    missing = [path for path in required_paths if not path.is_file()]
    if missing:
        raise RuntimeError(
            f"Cannot checkpoint {stage}; missing output(s): "
            + ", ".join(map(str, missing))
        )
    _atomic_json(
        marker,
        {
            "schema_version": 1,
            "complete": True,
            "stage": stage,
            "signature": signature,
            "sizes": {path.name: path.stat().st_size for path in required_paths},
        },
    )


def _stream_percentile_group(
    *,
    source_by_chrom: dict[str, list[PositionRecord]],
    target_by_chrom: dict[str, list[PositionRecord]],
    all_target_count: int,
    assignment_a: np.ndarray,
    assignment_b: np.ndarray,
    percentile_source: str,
    matching: str,
    max_distance: float | None,
    lower: int,
    upper: int,
    label: str,
    detail_path: Path | None,
    distance_path: Path,
    reporter: ProgressReporter,
) -> dict[str, object]:
    """Process one directional percentile group without retaining its pair rows."""

    target_source = "B" if percentile_source == "A" else "A"
    source_assignment = assignment_a if percentile_source == "A" else assignment_b
    detail_handle = None
    detail_writer = None
    temporary_detail = None
    if detail_path is not None:
        detail_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_detail = detail_path.with_name(detail_path.name + ".partial")
        detail_handle = temporary_detail.open("w", encoding="utf-8", newline="")
        detail_writer = csv.DictWriter(
            detail_handle, fieldnames=PERCENTILE_DETAIL_FIELDS, delimiter="\t"
        )
        detail_writer.writeheader()

    distances = array("q")
    source_count = 0
    matched_pairs = 0
    matched_unique_targets = 0
    unmatched_no_target = 0
    unmatched_distance = 0
    unmatched_unique = 0
    pair_number = 0
    contigs = list(source_by_chrom)
    try:
        for contig_index, (chrom_key, source_records) in enumerate(
            source_by_chrom.items(), start=1
        ):
            query_group = [
                record
                for record in source_records
                if lower < int(source_assignment[record.line_number]) <= upper
            ]
            if not query_group:
                continue
            source_count += len(query_group)
            reporter.contig(
                f"{percentile_source} {label} percentile matching",
                query_group[0].chrom,
                contig_index,
                len(contigs),
                len(query_group),
            )
            targets = target_by_chrom.get(chrom_key)
            if not targets:
                unmatched_no_target += len(query_group)
                continue
            if matching == "many-to-one":
                result = match_many_to_one(
                    query_group,
                    targets,
                    percentile_source,
                    max_distance,
                )
            else:
                result = match_unique(
                    query_group,
                    targets,
                    percentile_source,
                    max_distance,
                )
            matched_pairs += len(result.pairs)
            unmatched_no_target += result.unmatched_no_target_chrom
            unmatched_distance += result.unmatched_distance
            unmatched_unique += result.unmatched_unique
            matched_target_lines: set[int] = set()
            for pair in result.pairs:
                pair_number += 1
                distances.append(pair.absolute_distance)
                target_record = pair.b if percentile_source == "A" else pair.a
                matched_target_lines.add(target_record.line_number)
                if detail_writer is None:
                    continue
                a_percentile = int(assignment_a[pair.a.line_number])
                b_percentile = int(assignment_b[pair.b.line_number])
                signed_target_minus_query = (
                    pair.signed_distance
                    if percentile_source == "A"
                    else -pair.signed_distance
                )
                row = {
                    "pair_id": (
                        f"{percentile_source}_{label.replace('-', '_')}_"
                        f"{pair_number:09d}"
                    ),
                    "analysis_direction": (
                        f"{percentile_source}_percentiles_vs_all_{target_source}"
                    ),
                    "percentile_source": percentile_source,
                    "target_source": target_source,
                    "percentile_group": label,
                    "group_lower_percentile": lower,
                    "group_upper_percentile": upper,
                    "query_source": percentile_source,
                    "chrom": pair.a.chrom,
                    "a_start": pair.a.start,
                    "a_end": pair.a.end,
                    "a_name": pair.a.name,
                    "a_summit": pair.a.summit,
                    "a_score": pair.a.score,
                    "a_score_percentile": a_percentile,
                    "b_start": pair.b.start,
                    "b_end": pair.b.end,
                    "b_name": pair.b.name,
                    "b_summit": pair.b.summit,
                    "b_score": pair.b.score,
                    "b_score_percentile": b_percentile,
                    "signed_distance_b_minus_a": pair.signed_distance,
                    "signed_distance_target_minus_query": signed_target_minus_query,
                    "absolute_distance": pair.absolute_distance,
                }
                detail_writer.writerow(
                    {key: _format_number(row[key]) for key in PERCENTILE_DETAIL_FIELDS}
                )
            matched_unique_targets += len(matched_target_lines)
            del result
        distance_values = np.frombuffer(distances, dtype=np.int64).copy()
        _save_numpy_atomic(distance_path, distance_values)
        if distance_values.size:
            q1, median, q3 = np.percentile(distance_values, [25, 50, 75])
            minimum = float(np.min(distance_values))
            maximum = float(np.max(distance_values))
            mean = float(np.mean(distance_values))
            standard_deviation = float(np.std(distance_values, ddof=0))
        else:
            minimum = q1 = median = mean = q3 = maximum = standard_deviation = float("nan")
        summary = {
            "analysis_direction": f"{percentile_source}_percentiles_vs_all_{target_source}",
            "percentile_source": percentile_source,
            "target_source": target_source,
            "percentile_group": label,
            "group_lower_percentile": lower,
            "group_upper_percentile": upper,
            "percentile_source_position_count": source_count,
            "all_target_position_count": all_target_count,
            "matched_pair_count": matched_pairs,
            "matched_unique_source_positions": matched_pairs,
            "matched_unique_target_positions": matched_unique_targets,
            "unmatched_source_positions": source_count - matched_pairs,
            "target_positions_not_used": all_target_count - matched_unique_targets,
            "unmatched_query_no_target_chromosome": unmatched_no_target,
            "unmatched_query_beyond_maximum_distance": unmatched_distance,
            "unmatched_query_unique_assignment": unmatched_unique,
            "minimum_absolute_distance": minimum,
            "q1_absolute_distance": q1,
            "median_absolute_distance": median,
            "mean_absolute_distance": mean,
            "q3_absolute_distance": q3,
            "maximum_absolute_distance": maximum,
            "absolute_distance_standard_deviation": standard_deviation,
        }
    finally:
        if detail_handle is not None:
            detail_handle.close()
    if detail_path is not None and temporary_detail is not None:
        os.replace(temporary_detail, detail_path)
    return summary


def _concatenate_percentile_groups(inputs: Sequence[Path], output: Path) -> None:
    temporary = output.with_name(output.name + ".partial")
    output.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open("w", encoding="utf-8", newline="") as destination:
        wrote_header = False
        for path in inputs:
            with path.open("r", encoding="utf-8") as source:
                header = source.readline()
                if not wrote_header:
                    destination.write(header)
                    wrote_header = True
                for line in source:
                    destination.write(line)
    os.replace(temporary, output)


def run_directional_percentile_streaming(
    *,
    records_a: Sequence[PositionRecord],
    records_b: Sequence[PositionRecord],
    assignment_a: np.ndarray,
    assignment_b: np.ndarray,
    percentile_source: str,
    matching: str,
    max_distance: float | None,
    interval: int,
    output_distances: Path | None,
    output_summary: Path,
    output_plot: Path,
    checkpoint_root: Path,
    signature: str,
    label_source: str,
    label_target: str,
    dpi: int,
    y_max: float,
    force: bool,
    reporter: ProgressReporter,
) -> dict[str, Path]:
    source_records = records_a if percentile_source == "A" else records_b
    target_records = records_b if percentile_source == "A" else records_a
    source_by_chrom = _records_by_chrom(source_records)
    target_by_chrom = _records_by_chrom(target_records)
    group_details: list[Path] = []
    summary_rows: list[dict[str, object]] = []
    distances_by_group: dict[str, np.ndarray] = {}
    direction_dir = checkpoint_root / f"{percentile_source}_percentiles"
    direction_dir.mkdir(parents=True, exist_ok=True)

    for lower, upper, label in _percentile_group_bounds(interval):
        safe_label = label.replace("-", "_")
        detail_path = (
            direction_dir / f"{safe_label}.pairs.tsv"
            if output_distances is not None
            else None
        )
        distance_path = direction_dir / f"{safe_label}.distances.npy"
        summary_path = direction_dir / f"{safe_label}.summary.json"
        marker = direction_dir / f"{safe_label}.complete.json"
        group_signature = __import__("hashlib").sha256(
            f"{signature}|{percentile_source}|{lower}|{upper}|{output_distances is not None}".encode()
        ).hexdigest()
        required = [distance_path, summary_path]
        if detail_path is not None:
            required.append(detail_path)
        if not force and _checkpoint_matches(
            marker, signature=group_signature, required_paths=required
        ):
            reporter.emit(
                f"Reused completed {percentile_source} percentile group {label}"
            )
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        else:
            marker.unlink(missing_ok=True)
            reporter.stage(
                f"Processing {percentile_source} percentile group {label}"
            )
            summary = _stream_percentile_group(
                source_by_chrom=source_by_chrom,
                target_by_chrom=target_by_chrom,
                all_target_count=len(target_records),
                assignment_a=assignment_a,
                assignment_b=assignment_b,
                percentile_source=percentile_source,
                matching=matching,
                max_distance=max_distance,
                lower=lower,
                upper=upper,
                label=label,
                detail_path=detail_path,
                distance_path=distance_path,
                reporter=reporter,
            )
            _atomic_json(summary_path, summary)
            _atomic_json(
                marker,
                {
                    "schema_version": 1,
                    "complete": True,
                    "signature": group_signature,
                    "sizes": {
                        path.name: path.stat().st_size for path in required
                    },
                },
            )
            reporter.emit(
                f"Completed {percentile_source} percentile group {label}; "
                f"{int(summary['matched_pair_count']):,} pairs"
            )
        summary_rows.append(summary)
        distances_by_group[label] = np.load(
            distance_path, mmap_mode="r", allow_pickle=False
        )
        if detail_path is not None:
            group_details.append(detail_path)

    if output_distances is not None:
        _concatenate_percentile_groups(group_details, output_distances)
    write_percentile_distance_summary(output_summary, summary_rows)
    plot_percentile_distances(
        output_plot,
        distances_by_group,
        interval,
        label_source,
        label_target,
        dpi,
        y_max,
    )
    outputs = {"summary": output_summary, "boxplot": output_plot}
    if output_distances is not None:
        outputs["distances"] = output_distances
    return outputs


def _safe_zscores(values: np.ndarray) -> np.ndarray:
    if values.size == 0:
        return values.copy()
    standard_deviation = float(np.std(values, ddof=0))
    if standard_deviation == 0 or not math.isfinite(standard_deviation):
        return np.zeros(values.size, dtype=np.float64)
    return (values - float(np.mean(values))) / standard_deviation


def _percentile_ranks(values: np.ndarray) -> np.ndarray:
    if values.size == 0:
        return values.copy()
    if values.size == 1:
        return np.asarray([0.5], dtype=np.float64)
    ranks = stats.rankdata(values, method="average")
    return (ranks - 1.0) / (values.size - 1.0)


def _safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(stats.pearsonr(x, y).statistic)


def _safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(stats.spearmanr(x, y).statistic)


def _safe_regression(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    if x.size < 2 or np.ptp(x) == 0:
        return float("nan"), float("nan"), float("nan")
    regression = stats.linregress(x, y)
    return float(regression.slope), float(regression.intercept), float(regression.rvalue**2)


def _parse_distance_bins(text: str) -> list[float]:
    bounds: list[float] = []
    for token in text.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        try:
            value = float(stripped)
        except ValueError as exc:
            raise ValueError(f"Invalid distance-bin boundary: {stripped!r}.") from exc
        if not math.isfinite(value) or value < 0:
            raise ValueError("Distance-bin boundaries must be finite and non-negative.")
        bounds.append(value)
    bounds = sorted(set(bounds))
    if not bounds:
        raise ValueError("--distance-bins must contain at least one boundary.")
    return bounds


def _bin_label(lower: float | None, upper: float | None) -> str:
    if lower is None:
        return f"0-{upper:g}"
    if upper is None:
        return f">{lower:g}"
    if float(lower).is_integer() and float(upper).is_integer():
        return f"{int(lower) + 1}-{int(upper)}"
    return f"{lower:g}-{upper:g}"


def _selected_scores(arrays: dict[str, np.ndarray], normalization: str) -> tuple[np.ndarray, np.ndarray, str]:
    if normalization == "raw":
        return arrays["scores_a"], arrays["scores_b"], "raw score"
    if normalization == "zscore":
        return arrays["z_a"], arrays["z_b"], "score z-score"
    if normalization == "percentile":
        return arrays["percentile_a"], arrays["percentile_b"], "score percentile rank"
    raise ValueError(f"Unknown score normalization: {normalization}")


def distance_bin_rows(
    distances: np.ndarray,
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    z_difference: np.ndarray,
    raw_difference: np.ndarray,
    bounds: Sequence[float],
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    lower: float | None = None
    for upper in [*bounds, None]:
        if lower is None:
            mask = distances <= float(upper)
        elif upper is None:
            mask = distances > lower
        else:
            mask = (distances > lower) & (distances <= upper)
        bin_a = scores_a[mask]
        bin_b = scores_b[mask]
        bin_distances = distances[mask]
        rows.append(
            {
                "distance_bin": _bin_label(lower, upper),
                "lower_exclusive": "" if lower is None else lower,
                "upper_inclusive": "" if upper is None else upper,
                "pair_count": int(np.sum(mask)),
                "mean_absolute_distance": float(np.mean(bin_distances)) if bin_distances.size else float("nan"),
                "median_absolute_distance": float(np.median(bin_distances)) if bin_distances.size else float("nan"),
                "pearson_score_correlation": _safe_pearson(bin_a, bin_b),
                "spearman_score_correlation": _safe_spearman(bin_a, bin_b),
                "median_absolute_score_difference": float(np.median(raw_difference[mask])) if bin_distances.size else float("nan"),
                "median_absolute_zscore_difference": float(np.median(z_difference[mask])) if bin_distances.size else float("nan"),
            }
        )
        lower = upper
    return rows


def _format_number(value: object) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        return f"{value:.10g}"
    return str(value)


def write_pairs(
    path: Path | None,
    pairs: Sequence[MatchedPair],
    *,
    plot_indices: np.ndarray | None = None,
    plot_score_normalization: str = "zscore",
    plot_correlation_method: str = "spearman",
    plot_label_a: str = "A",
    plot_label_b: str = "B",
    plot_score_z_limit: float = 10.0,
) -> dict[str, np.ndarray]:
    """Create analysis arrays and optionally publish the detailed pair table."""
    scores_a = np.asarray([pair.a.score for pair in pairs], dtype=np.float64)
    scores_b = np.asarray([pair.b.score for pair in pairs], dtype=np.float64)
    distances = np.asarray([pair.absolute_distance for pair in pairs], dtype=np.float64)
    z_a = _safe_zscores(scores_a)
    z_b = _safe_zscores(scores_b)
    percentile_a = _percentile_ranks(scores_a)
    percentile_b = _percentile_ranks(scores_b)
    raw_difference = np.abs(scores_b - scores_a)
    z_difference = np.abs(z_b - z_a)
    percentile_difference = np.abs(percentile_b - percentile_a)

    fieldnames = [
        "pair_id",
        "query_method",
        "chrom",
        "a_start",
        "a_end",
        "a_name",
        "a_summit",
        "a_score",
        "b_start",
        "b_end",
        "b_name",
        "b_summit",
        "b_score",
        "signed_distance_b_minus_a",
        "absolute_distance",
        "score_difference_b_minus_a",
        "absolute_score_difference",
        "a_score_z",
        "b_score_z",
        "absolute_zscore_difference",
        "a_score_percentile",
        "b_score_percentile",
        "absolute_percentile_difference",
        "plot_selected",
        "plot_score_normalization",
        "plot_correlation_method",
        "plot_label_a",
        "plot_label_b",
        "plot_score_z_limit",
    ]
    selected = (
        set(int(value) for value in plot_indices)
        if plot_indices is not None else set(range(len(pairs)))
    )
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".partial")
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for index, pair in enumerate(pairs):
                writer.writerow(
                    {
                        "pair_id": f"pair_{index + 1:07d}",
                        "query_method": pair.query_source,
                        "chrom": pair.a.chrom,
                        "a_start": pair.a.start,
                        "a_end": pair.a.end,
                        "a_name": pair.a.name,
                        "a_summit": pair.a.summit,
                        "a_score": _format_number(pair.a.score),
                        "b_start": pair.b.start,
                        "b_end": pair.b.end,
                        "b_name": pair.b.name,
                        "b_summit": pair.b.summit,
                        "b_score": _format_number(pair.b.score),
                        "signed_distance_b_minus_a": pair.signed_distance,
                        "absolute_distance": pair.absolute_distance,
                        "score_difference_b_minus_a": _format_number(pair.b.score - pair.a.score),
                        "absolute_score_difference": _format_number(raw_difference[index]),
                        "a_score_z": _format_number(z_a[index]),
                        "b_score_z": _format_number(z_b[index]),
                        "absolute_zscore_difference": _format_number(z_difference[index]),
                        "a_score_percentile": _format_number(percentile_a[index]),
                        "b_score_percentile": _format_number(percentile_b[index]),
                        "absolute_percentile_difference": _format_number(percentile_difference[index]),
                        "plot_selected": 1 if index in selected else 0,
                        "plot_score_normalization": plot_score_normalization,
                        "plot_correlation_method": plot_correlation_method,
                        "plot_label_a": plot_label_a,
                        "plot_label_b": plot_label_b,
                        "plot_score_z_limit": _format_number(plot_score_z_limit),
                    }
                )
        os.replace(temporary, path)
    return {
        "scores_a": scores_a,
        "scores_b": scores_b,
        "distances": distances,
        "z_a": z_a,
        "z_b": z_b,
        "percentile_a": percentile_a,
        "percentile_b": percentile_b,
        "raw_difference": raw_difference,
        "z_difference": z_difference,
        "percentile_difference": percentile_difference,
    }


def summary_metrics(
    result: MatchResult,
    arrays: dict[str, np.ndarray],
    matching: str,
    max_distance: float | None,
    records_a: Sequence[PositionRecord],
    records_b: Sequence[PositionRecord],
    score_normalization: str,
    correlation_method: str,
) -> list[tuple[str, object]]:
    scores_a, scores_b, _ = _selected_scores(arrays, score_normalization)
    distances = arrays["distances"]
    slope, intercept, r_squared = _safe_regression(scores_a, scores_b)
    metrics: list[tuple[str, object]] = [
        ("input_a_positions", len(records_a)),
        ("input_b_positions", len(records_b)),
        ("query_method", result.query_source),
        ("target_method", result.target_source),
        ("matching_policy", matching),
        ("maximum_distance_bp", "none" if max_distance is None else max_distance),
        ("score_normalization", score_normalization),
        ("correlation_method", correlation_method),
        ("matched_pairs", len(result.pairs)),
        ("unmatched_no_target_chromosome", result.unmatched_no_target_chrom),
        ("unmatched_beyond_maximum_distance", result.unmatched_distance),
        ("unmatched_unique_assignment", result.unmatched_unique),
    ]
    if distances.size:
        metrics.extend(
            [
                ("mean_signed_distance_b_minus_a", float(np.mean([pair.signed_distance for pair in result.pairs]))),
                ("median_signed_distance_b_minus_a", float(np.median([pair.signed_distance for pair in result.pairs]))),
                ("mean_absolute_distance", float(np.mean(distances))),
                ("median_absolute_distance", float(np.median(distances))),
                ("absolute_distance_standard_deviation", float(np.std(distances))),
                ("summit_distance_rmse", float(np.sqrt(np.mean(np.square([pair.signed_distance for pair in result.pairs]))))),
            ]
        )
        if correlation_method in ("spearman", "both"):
            metrics.append(("spearman_score_correlation", _safe_spearman(scores_a, scores_b)))
        if correlation_method in ("pearson", "both"):
            metrics.append(("pearson_score_correlation", _safe_pearson(scores_a, scores_b)))
        metrics.extend(
            [
                ("score_regression_slope_b_on_a", slope),
                ("score_regression_intercept_b_on_a", intercept),
                ("score_regression_r_squared", r_squared),
            ]
        )
        for threshold in (5, 10, 20, 50, 100):
            metrics.append((f"percent_pairs_within_{threshold}_bp", float(100.0 * np.mean(distances <= threshold))))
    return metrics


def write_summary(path: Path, metrics: Sequence[tuple[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for metric, value in metrics:
            handle.write(f"{metric}\t{_format_number(value)}\n")
    os.replace(temporary, path)


def _update_summary_metric(path: Path, metric: str, value: object) -> None:
    """Atomically replace or append one metric in a two-column summary."""
    rows: list[tuple[str, str]] = []
    found = False
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if header.rstrip("\n") != "metric\tvalue":
            raise ValueError(f"Unexpected summary header in {path}")
        for raw in handle:
            key, old_value = raw.rstrip("\n").split("\t", 1)
            if key == metric:
                rows.append((key, _format_number(value)))
                found = True
            else:
                rows.append((key, old_value))
    if not found:
        rows.append((metric, _format_number(value)))
    write_summary(path, rows)


def write_distance_bins(
    path: Path,
    rows: Sequence[dict[str, object]],
    *,
    plot_correlation_method: str = "spearman",
    plot_score_axis_label: str = "score z-score",
) -> None:
    fieldnames = [
        "distance_bin",
        "lower_exclusive",
        "upper_inclusive",
        "pair_count",
        "mean_absolute_distance",
        "median_absolute_distance",
        "pearson_score_correlation",
        "spearman_score_correlation",
        "median_absolute_score_difference",
        "median_absolute_zscore_difference",
        "plot_correlation_method",
        "plot_score_axis_label",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            output_row = {key: _format_number(row[key]) for key in fieldnames if key in row}
            output_row["plot_correlation_method"] = plot_correlation_method
            output_row["plot_score_axis_label"] = plot_score_axis_label
            writer.writerow(output_row)



def _histogram_edges(x_min: float, x_max: float, bin_width: float) -> np.ndarray:
    """Return histogram edges spanning the requested displayed range."""

    edges = np.arange(x_min, x_max + bin_width, bin_width, dtype=np.float64)
    if edges.size < 2:
        edges = np.asarray([x_min, x_max], dtype=np.float64)
    elif edges[-1] < x_max:
        edges = np.append(edges, x_max)
    elif edges[-1] > x_max:
        edges[-1] = x_max
    return edges


def write_distance_histogram(
    path: Path,
    distances: np.ndarray,
    bin_width: float,
    x_min: float,
    x_max: float,
) -> None:
    """Write the exact binned counts used by the distance histogram plot."""

    edges = _histogram_edges(x_min, x_max, bin_width)
    counts, _ = np.histogram(distances, bins=edges)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "bin_start_inclusive",
            "bin_end_exclusive",
            "distance_bin",
            "pair_count",
        ])
        for index, count in enumerate(counts):
            start = float(edges[index])
            end = float(edges[index + 1])
            writer.writerow([
                _format_number(start),
                _format_number(end),
                f"{_format_number(start)}-{_format_number(end)}",
                int(count),
            ])

def write_percentile_distances(path: Path, rows: Sequence[dict[str, object]]) -> None:
    fieldnames = [
        "pair_id",
        "analysis_direction",
        "percentile_source",
        "target_source",
        "percentile_group",
        "group_lower_percentile",
        "group_upper_percentile",
        "query_source",
        "chrom",
        "a_start",
        "a_end",
        "a_name",
        "a_summit",
        "a_score",
        "a_score_percentile",
        "b_start",
        "b_end",
        "b_name",
        "b_summit",
        "b_score",
        "b_score_percentile",
        "signed_distance_b_minus_a",
        "signed_distance_target_minus_query",
        "absolute_distance",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format_number(row[key]) for key in fieldnames})


def write_percentile_distance_summary(
    path: Path,
    rows: Sequence[dict[str, object]],
) -> None:
    fieldnames = [
        "analysis_direction",
        "percentile_source",
        "target_source",
        "percentile_group",
        "group_lower_percentile",
        "group_upper_percentile",
        "percentile_source_position_count",
        "all_target_position_count",
        "matched_pair_count",
        "matched_unique_source_positions",
        "matched_unique_target_positions",
        "unmatched_source_positions",
        "target_positions_not_used",
        "unmatched_query_no_target_chromosome",
        "unmatched_query_beyond_maximum_distance",
        "unmatched_query_unique_assignment",
        "minimum_absolute_distance",
        "q1_absolute_distance",
        "median_absolute_distance",
        "mean_absolute_distance",
        "q3_absolute_distance",
        "maximum_absolute_distance",
        "absolute_distance_standard_deviation",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format_number(row[key]) for key in fieldnames})


def plot_percentile_distances(
    path: Path,
    distances_by_group: dict[str, np.ndarray],
    interval: int,
    percentile_label: str,
    target_label: str,
    dpi: int,
    y_max: float = 500.0,
) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    groups = _percentile_group_bounds(interval)
    labels = [label for _, _, label in groups]
    figure_width = max(8.5, 0.72 * len(labels) + 2.5)
    figure, axis = plt.subplots(figsize=(figure_width, 6.2))
    for position, label in enumerate(labels, start=1):
        values = distances_by_group.get(label, np.asarray([], dtype=np.float64))
        if values.size:
            axis.boxplot(
                [values],
                positions=[position],
                widths=0.58,
                showfliers=True,
                whis=1.5,
            )
    axis.set_xticks(np.arange(1, len(labels) + 1), labels, rotation=35, ha="right")
    axis.set_xlim(0.4, len(labels) + 0.6)
    axis.set_xlabel(f"{percentile_label} score percentile group")
    axis.set_ylabel("Absolute summit distance (bp)")
    if y_max > 0:
        axis.set_ylim(0.0, y_max)
    axis.set_title(
        f"{percentile_label} score percentiles versus all {target_label} positions"
    )
    from nucleosuite.plotting import save_figure
    figure.tight_layout()
    saved = save_figure(figure, path, default_dpi=dpi)
    plt.close(figure)
    return saved


def _plot_indices(pair_count: int, maximum: int, seed: int) -> np.ndarray:
    if maximum == 0 or pair_count <= maximum:
        return np.arange(pair_count, dtype=int)
    random_generator = np.random.default_rng(seed)
    return np.sort(random_generator.choice(pair_count, size=maximum, replace=False))


def create_plots(
    prefix: Path,
    arrays: dict[str, np.ndarray],
    bin_rows: Sequence[dict[str, object]],
    label_a: str,
    label_b: str,
    plot_max_points: int,
    plot_seed: int,
    dpi: int,
    score_normalization: str,
    correlation_method: str,
    histogram_bin_width: float,
    histogram_x_min: float,
    histogram_x_max: float,
    score_z_limit: float = 10.0,
    distance_x_major_tick: float | None = None,
    distance_x_minor_tick: float | None = None,
) -> list[Path]:
    if plot_max_points < 0:
        raise ValueError("--plot-max-points must be zero or greater.")
    if dpi < 1:
        raise ValueError("--dpi must be at least 1.")
    if histogram_bin_width <= 0:
        raise ValueError("--histogram-bin-width must be greater than zero.")
    if histogram_x_min < 0 or histogram_x_max <= histogram_x_min:
        raise ValueError("Histogram limits require 0 <= x-min < x-max.")
    if not math.isfinite(score_z_limit) or score_z_limit < 0:
        raise ValueError("--score-z-limit must be finite and zero or greater.")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    scores_a, scores_b, score_axis_label = _selected_scores(arrays, score_normalization)
    distances = arrays["distances"]
    indices = _plot_indices(scores_a.size, plot_max_points, plot_seed)
    outputs: list[Path] = []

    from nucleosuite.plotting import plot_path, save_figure
    score_plot = plot_path(Path(f"{prefix}_score_correlation.png"))
    figure, axis = plt.subplots(figsize=(7.5, 6.5))
    scatter = axis.scatter(scores_a[indices], scores_b[indices], c=distances[indices], s=9, alpha=0.55, linewidths=0, cmap="viridis", rasterized=True)
    colorbar = figure.colorbar(scatter, ax=axis)
    colorbar.set_label("Absolute summit distance (bp)")
    _slope, _intercept, r_squared = _safe_regression(scores_a, scores_b)
    if score_normalization == "zscore" and score_z_limit > 0:
        axis.set_xlim(-score_z_limit, score_z_limit)
        axis.set_ylim(-score_z_limit, score_z_limit)
    axis.set_xlabel(f"{label_a} {score_axis_label}")
    axis.set_ylabel(f"{label_b} {score_axis_label}")
    axis.set_title("Score agreement coloured by summit distance")
    annotation = []
    if correlation_method in ("spearman", "both"):
        annotation.append(f"Spearman ρ = {_safe_spearman(scores_a, scores_b):.3f}")
    if correlation_method in ("pearson", "both"):
        annotation.append(f"Pearson r = {_safe_pearson(scores_a, scores_b):.3f}")
    annotation.append(f"R² = {r_squared:.3f}")
    axis.text(0.02, 0.98, "\n".join(annotation), transform=axis.transAxes, va="top", ha="left")
    figure.tight_layout()
    score_plot = save_figure(figure, score_plot, default_dpi=dpi)
    plt.close(figure)
    outputs.append(score_plot)

    histogram_plot = plot_path(Path(f"{prefix}_distance_histogram.png"))
    figure, axis = plt.subplots(figsize=(7.5, 5.5))
    bin_edges = _histogram_edges(histogram_x_min, histogram_x_max, histogram_bin_width)
    axis.hist(distances, bins=bin_edges)
    for threshold in (5, 10, 20, 50, 100):
        if histogram_x_min <= threshold <= histogram_x_max:
            axis.axvline(threshold, color="black", linewidth=0.7, linestyle="--", alpha=0.45)
    axis.set_xlim(histogram_x_min, histogram_x_max)
    axis.set_xlabel("Absolute summit distance (bp)")
    axis.set_ylabel("Matched pairs")
    axis.set_title("Nearest-summit distance distribution")
    from nucleosuite.plotting import apply_distance_x_axis
    apply_distance_x_axis(
        axis,
        major_interval=distance_x_major_tick,
        minor_interval=distance_x_minor_tick,
    )
    from nucleosuite.plotting import apply_integer_y_axis
    apply_integer_y_axis(axis)
    figure.tight_layout()
    histogram_plot = save_figure(figure, histogram_plot, default_dpi=dpi)
    plt.close(figure)
    outputs.append(histogram_plot)

    bin_plot = plot_path(Path(f"{prefix}_correlation_by_distance.png"))
    figure, axis = plt.subplots(figsize=(8.5, 5.8))
    labels = [str(row["distance_bin"]) for row in bin_rows]
    x_values = np.arange(len(labels))
    if correlation_method in ("spearman", "both"):
        values = np.asarray([float(row["spearman_score_correlation"]) for row in bin_rows], dtype=np.float64)
        axis.plot(x_values, values, marker="o", label="Spearman")
    if correlation_method in ("pearson", "both"):
        values = np.asarray([float(row["pearson_score_correlation"]) for row in bin_rows], dtype=np.float64)
        axis.plot(x_values, values, marker="s", label="Pearson")
    axis.axhline(0, color="black", linewidth=0.8)
    axis.set_xticks(x_values, labels, rotation=35, ha="right")
    axis.set_ylim(-1.05, 1.05)
    axis.set_xlabel("Absolute summit-distance bin (bp)")
    axis.set_ylabel(f"{score_axis_label.capitalize()} correlation")
    axis.set_title("Score correlation by summit-distance bin")
    if correlation_method == "both":
        axis.legend()
    for x_value, row in zip(x_values, bin_rows):
        axis.text(x_value, -1.0, f"n={int(row['pair_count'])}", ha="center", va="bottom", fontsize=8, rotation=90)
    figure.tight_layout()
    bin_plot = save_figure(figure, bin_plot, default_dpi=dpi)
    plt.close(figure)
    outputs.append(bin_plot)
    return outputs


def _default_label(path: str | Path) -> str:
    name = Path(path).name
    for suffix in (".bed.gz", ".bigBed", ".bigbed", ".bed", ".bb", ".gz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _default_prefix(path_a: str | Path, path_b: str | Path) -> Path:
    return Path(f"{_default_label(path_a)}_vs_{_default_label(path_b)}")


def run_comparison(args: argparse.Namespace) -> dict[str, Path]:
    if args.plot_max_points < 0:
        raise ValueError("--plot-max-points must be zero or greater.")
    from nucleosuite.plotting import validate_tick_interval
    distance_x_major_tick = getattr(args, "distance_x_major_tick", None)
    distance_x_minor_tick = getattr(args, "distance_x_minor_tick", None)
    validate_tick_interval(distance_x_major_tick, "--distance-x-major-tick")
    validate_tick_interval(distance_x_minor_tick, "--distance-x-minor-tick")
    score_z_limit = float(getattr(args, "score_z_limit", 10.0))
    percentile_boxplot_y_max = float(
        getattr(args, "percentile_boxplot_y_max", 500.0)
    )
    if not math.isfinite(score_z_limit) or score_z_limit < 0:
        raise ValueError("--score-z-limit must be finite and zero or greater.")
    if not math.isfinite(percentile_boxplot_y_max) or percentile_boxplot_y_max < 0:
        raise ValueError(
            "--percentile-boxplot-y-max must be finite and zero or greater."
        )
    percentile_interval = _validate_percentile_interval(
        int(getattr(args, "percentile_interval", 25))
    )
    skip_percentile_analysis = bool(
        getattr(args, "skip_percentile_distance_analysis", False)
    )
    prefix = Path(args.output_prefix) if args.output_prefix else _default_prefix(args.bed_a, args.bed_b)
    from nucleosuite.output_naming import parameterized_prefix

    prefix = parameterized_prefix(
        prefix,
        (
            ("match", "one-to-one" if args.matching == "unique" else args.matching),
            ("maxdist", args.max_distance),
            ("scorenorm", args.score_normalization),
        ),
    )
    prefix.parent.mkdir(parents=True, exist_ok=True)
    reporter = ProgressReporter(
        "compare-positions", quiet=bool(getattr(args, "quiet", False))
    )
    bounds = _parse_distance_bins(args.distance_bins)
    pairs_path = Path(f"{prefix}_pairs.tsv")
    summary_path = Path(f"{prefix}_summary.tsv")
    bins_path = Path(f"{prefix}_distance_bins.tsv")
    histogram_tsv_path = Path(f"{prefix}_distance_histogram.tsv")
    a_percentile_distances_path = Path(
        f"{prefix}_A_percentiles_vs_all_B_distances.tsv"
    )
    a_percentile_summary_path = Path(
        f"{prefix}_A_percentiles_vs_all_B_summary.tsv"
    )
    from nucleosuite.plotting import plot_path
    a_percentile_plot_path = plot_path(Path(
        f"{prefix}_A_percentiles_vs_all_B_boxplot.png"
    ))
    b_percentile_distances_path = Path(
        f"{prefix}_B_percentiles_vs_all_A_distances.tsv"
    )
    b_percentile_summary_path = Path(
        f"{prefix}_B_percentiles_vs_all_A_summary.tsv"
    )
    b_percentile_plot_path = plot_path(Path(
        f"{prefix}_B_percentiles_vs_all_A_boxplot.png"
    ))
    score_plot_path = plot_path(Path(f"{prefix}_score_correlation.png"))
    distance_plot_path = plot_path(Path(f"{prefix}_distance_histogram.png"))
    correlation_plot_path = plot_path(Path(f"{prefix}_correlation_by_distance.png"))
    checkpoint_root = Path(f"{prefix}_checkpoints")
    checkpoint_root.mkdir(parents=True, exist_ok=True)
    main_marker = checkpoint_root / "main.complete.json"
    signature = _comparison_signature(args)
    skip_pairs = bool(getattr(args, "skip_pairs_tsv", False))
    skip_percentile_pairs = bool(
        getattr(args, "skip_percentile_pairs_tsv", False)
    )
    force = bool(getattr(args, "force", False))

    blacklist = load_blacklist_unbounded(getattr(args, "blacklist_bed", None))
    excluded_a = [0]
    excluded_b = [0]
    records_a = read_positions(
        args.bed_a,
        "A",
        args.summit_column_a,
        args.score_column_a,
        blacklist,
        excluded_a,
        progress=reporter,
    )
    records_b = read_positions(
        args.bed_b,
        "B",
        args.summit_column_b,
        args.score_column_b,
        blacklist,
        excluded_b,
        progress=reporter,
    )

    main_required = [
        summary_path,
        bins_path,
        histogram_tsv_path,
        score_plot_path,
        distance_plot_path,
        correlation_plot_path,
    ]
    if not skip_pairs:
        main_required.append(pairs_path)

    label_a = args.label_a or _default_label(args.bed_a)
    label_b = args.label_b or _default_label(args.bed_b)
    if not force and _checkpoint_matches(
        main_marker, signature=signature, required_paths=main_required
    ):
        reporter.emit("Reused completed main comparison outputs")
    else:
        main_marker.unlink(missing_ok=True)
        reporter.stage(
            f"Matching {args.matching}: A={len(records_a):,}; B={len(records_b):,}"
        )
        result = match_positions(
            records_a,
            records_b,
            args.matching,
            args.max_distance,
            progress=reporter,
        )
        if not result.pairs:
            raise ValueError("No matched position pairs satisfied the selected criteria.")
        reporter.stage(
            f"Preparing main analyses from {len(result.pairs):,} matched pairs"
        )
        arrays = write_pairs(
            None if skip_pairs else pairs_path,
            result.pairs,
            plot_indices=_plot_indices(len(result.pairs), args.plot_max_points, args.plot_seed),
            plot_score_normalization=args.score_normalization,
            plot_correlation_method=args.correlation_method,
            plot_label_a=label_a,
            plot_label_b=label_b,
            plot_score_z_limit=score_z_limit,
        )
        metrics = [
            ("input_a_file", str(Path(args.bed_a))),
            ("input_b_file", str(Path(args.bed_b))),
            ("blacklist_bed", getattr(args, "blacklist_bed", None) or ""),
            ("blacklist_overlapping_a_records_excluded", excluded_a[0]),
            ("blacklist_overlapping_b_records_excluded", excluded_b[0]),
            ("method_a_label", label_a),
            ("method_b_label", label_b),
            (
                "percentile_distance_analysis",
                "disabled" if skip_percentile_analysis else "in_progress",
            ),
            ("score_percentile_interval", percentile_interval),
            ("score_z_axis_limit", score_z_limit),
            ("distance_histogram_x_min", args.histogram_x_min),
            ("distance_histogram_x_max", args.histogram_x_max),
            ("percentile_boxplot_y_max", percentile_boxplot_y_max),
            ("pairs_tsv", "disabled" if skip_pairs else "written"),
            (
                "percentile_pairs_tsv",
                "disabled" if skip_percentile_pairs else "written",
            ),
        ] + summary_metrics(
            result,
            arrays,
            args.matching,
            args.max_distance,
            records_a,
            records_b,
            args.score_normalization,
            args.correlation_method,
        )
        write_summary(summary_path, metrics)
        selected_a, selected_b, _ = _selected_scores(
            arrays, args.score_normalization
        )
        rows = distance_bin_rows(
            arrays["distances"],
            selected_a,
            selected_b,
            arrays["z_difference"],
            arrays["raw_difference"],
            bounds,
        )
        write_distance_bins(
            bins_path,
            rows,
            plot_correlation_method=args.correlation_method,
            plot_score_axis_label=_selected_scores(arrays, args.score_normalization)[2],
        )
        write_distance_histogram(
            histogram_tsv_path,
            arrays["distances"],
            args.histogram_bin_width,
            args.histogram_x_min,
            args.histogram_x_max,
        )
        create_plots(
            prefix,
            arrays,
            rows,
            label_a,
            label_b,
            args.plot_max_points,
            args.plot_seed,
            args.dpi,
            args.score_normalization,
            args.correlation_method,
            args.histogram_bin_width,
            args.histogram_x_min,
            args.histogram_x_max,
            score_z_limit,
            distance_x_major_tick,
            distance_x_minor_tick,
        )
        _write_checkpoint(
            main_marker,
            signature=signature,
            required_paths=main_required,
            stage="main comparison",
        )
        reporter.emit("Completed main comparison outputs")
        del selected_a, selected_b, rows, arrays, result
        gc.collect()

    percentile_outputs: dict[str, Path] = {}
    if not skip_percentile_analysis:
        _update_summary_metric(
            summary_path, "percentile_distance_analysis", "in_progress"
        )
        _write_checkpoint(
            main_marker,
            signature=signature,
            required_paths=main_required,
            stage="main comparison",
        )
        reporter.stage("Assigning independent score percentiles")
        assignment_a = assign_score_percentiles_compact(
            records_a, percentile_interval
        )
        assignment_b = assign_score_percentiles_compact(
            records_b, percentile_interval
        )
        run_directional_percentile_streaming(
            records_a=records_a,
            records_b=records_b,
            assignment_a=assignment_a,
            assignment_b=assignment_b,
            percentile_source="A",
            matching=args.matching,
            max_distance=args.max_distance,
            interval=percentile_interval,
            output_distances=(
                None if skip_percentile_pairs else a_percentile_distances_path
            ),
            output_summary=a_percentile_summary_path,
            output_plot=a_percentile_plot_path,
            checkpoint_root=checkpoint_root,
            signature=signature,
            label_source=label_a,
            label_target=label_b,
            dpi=args.dpi,
            y_max=percentile_boxplot_y_max,
            force=force,
            reporter=reporter,
        )
        gc.collect()
        run_directional_percentile_streaming(
            records_a=records_a,
            records_b=records_b,
            assignment_a=assignment_a,
            assignment_b=assignment_b,
            percentile_source="B",
            matching=args.matching,
            max_distance=args.max_distance,
            interval=percentile_interval,
            output_distances=(
                None if skip_percentile_pairs else b_percentile_distances_path
            ),
            output_summary=b_percentile_summary_path,
            output_plot=b_percentile_plot_path,
            checkpoint_root=checkpoint_root,
            signature=signature,
            label_source=label_b,
            label_target=label_a,
            dpi=args.dpi,
            y_max=percentile_boxplot_y_max,
            force=force,
            reporter=reporter,
        )
        _update_summary_metric(
            summary_path, "percentile_distance_analysis", "completed"
        )
        _write_checkpoint(
            main_marker,
            signature=signature,
            required_paths=main_required,
            stage="main comparison",
        )
        percentile_outputs = {
            "a_percentiles_vs_all_b_summary": a_percentile_summary_path,
            "a_percentiles_vs_all_b_boxplot": a_percentile_plot_path,
            "b_percentiles_vs_all_a_summary": b_percentile_summary_path,
            "b_percentiles_vs_all_a_boxplot": b_percentile_plot_path,
        }
        if not skip_percentile_pairs:
            percentile_outputs.update(
                {
                    "a_percentiles_vs_all_b_distances": a_percentile_distances_path,
                    "b_percentiles_vs_all_a_distances": b_percentile_distances_path,
                }
            )

    outputs = {
        "summary": summary_path,
        "distance_bins": bins_path,
        "distance_histogram": histogram_tsv_path,
        "score_correlation_plot": score_plot_path,
        "distance_histogram_plot": distance_plot_path,
        "correlation_by_distance_plot": correlation_plot_path,
        **percentile_outputs,
    }
    if not skip_pairs:
        outputs["pairs"] = pairs_path
    if not getattr(args, "quiet", False):
        reporter.emit(
            f"A={len(records_a):,}; B={len(records_b):,}; "
            f"blacklisted A/B={excluded_a[0]:,}/{excluded_b[0]:,}"
        )
        for name, path in outputs.items():
            reporter.emit(f"{name}: {path}")
    reporter.complete()
    return outputs


def _run_serial(args: argparse.Namespace) -> int:
    try:
        run_comparison(args)
    except (FileNotFoundError, OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


def run(args: argparse.Namespace) -> int:
    from nucleosuite.partitioned import run_partitioned_command
    from nucleosuite.output_naming import parameterized_prefix

    requested = args.output_prefix or _default_prefix(args.bed_a, args.bed_b)
    args.output_prefix = str(
        parameterized_prefix(
            requested,
            (
                ("match", "one-to-one" if args.matching == "unique" else args.matching),
                ("maxdist", args.max_distance),
                ("scorenorm", args.score_normalization),
            ),
        )
    )
    prefix = Path(args.output_prefix).name if args.output_prefix else _default_prefix(args.bed_a, args.bed_b).name
    return run_partitioned_command(
        "compare-positions", args, _run_serial,
        runner_module="nucleosuite.compare_positions", runner_function="_run_serial",
        primary_attr="bed_a", output_prefix_attr="output_prefix",
        path_attrs=("bed_b", "blacklist_bed"), base_name=prefix,
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
