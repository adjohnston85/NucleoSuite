"""Fragment-coordinate randomisation utilities.

The materialised workflow divides each contig into fixed local blocks. Candidate
coordinates are indexed only inside the core block; BAM/reference flanks used to
recover complete source fragments are not candidate positions.
"""

from __future__ import annotations

import random
from array import array
from bisect import bisect_left
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

from nucleosuite.core.blacklist import BlacklistIndex

DINUCS = tuple(a + b for a in "ACGT" for b in "ACGT")
Fragment = Tuple[int, int]


def build_dinuc_index(sequence: str) -> Dict[str, List[int]]:
    sequence = sequence.upper()
    index = {dinuc: [] for dinuc in DINUCS}
    for position in range(len(sequence) - 1):
        dinuc = sequence[position : position + 2]
        if dinuc in index:
            index[dinuc].append(position)
    return index


def _invalid_prefix(sequence: str) -> list[int]:
    prefix = [0]
    total = 0
    for base in sequence.upper():
        total += base not in "ACGT"
        prefix.append(total)
    return prefix


@dataclass
class RandomizationBlock:
    """Candidate coordinates for one fixed half-open genomic block."""

    chrom: str
    start: int
    end: int
    sequence: str
    blacklist: BlacklistIndex | None = None
    coordinate_counts: MutableMapping[Fragment, int] | None = None
    max_per_coordinate: int = 0
    _dinuc_index: dict[str, list[int]] = field(init=False, repr=False)
    _invalid: list[int] = field(init=False, repr=False)
    _uniform_cache: dict[int, array] = field(default_factory=dict, init=False, repr=False)
    _anchor_cache: dict[tuple[str, int, str], array] = field(default_factory=dict, init=False, repr=False)
    blacklist_candidate_rejections: int = field(default=0, init=False)
    non_acgt_candidate_rejections: int = field(default=0, init=False)

    def __post_init__(self) -> None:
        if self.end <= self.start:
            raise ValueError("Randomization block must have positive length")
        if len(self.sequence) != self.end - self.start:
            raise ValueError("Randomization block sequence length does not match coordinates")
        self.sequence = self.sequence.upper()
        self._dinuc_index = build_dinuc_index(self.sequence)
        self._invalid = _invalid_prefix(self.sequence)
        if self.coordinate_counts is None:
            self.coordinate_counts = Counter()
        if self.max_per_coordinate < 0:
            raise ValueError("max_per_coordinate must be 0 or greater")

    @property
    def length(self) -> int:
        return self.end - self.start

    def _sequence_is_canonical(self, relative_start: int, relative_end: int) -> bool:
        return self._invalid[relative_end] == self._invalid[relative_start]

    def _start_is_valid(self, candidate_start: int, length: int) -> bool:
        candidate_end = candidate_start + length
        if candidate_start < self.start or candidate_end > self.end:
            return False
        relative_start = candidate_start - self.start
        if not self._sequence_is_canonical(relative_start, relative_start + length):
            self.non_acgt_candidate_rejections += 1
            return False
        if self.blacklist is not None and self.blacklist.overlaps(
            self.chrom, candidate_start, candidate_end
        ):
            self.blacklist_candidate_rejections += 1
            return False
        return True

    def _coordinate_available(self, fragment: Fragment) -> bool:
        if self.max_per_coordinate == 0:
            return True
        assert self.coordinate_counts is not None
        return self.coordinate_counts.get(fragment, 0) < self.max_per_coordinate

    def uniform_candidates(self, length: int, original: Fragment) -> list[int]:
        cached = self._uniform_starts(length)
        output: list[int] = []
        for candidate_start in cached:
            candidate = (int(candidate_start), int(candidate_start) + length)
            if candidate != original and self._coordinate_available(candidate):
                output.append(int(candidate_start))
        return output

    def _uniform_starts(self, length: int) -> array:
        cached = self._uniform_cache.get(length)
        if cached is None:
            maximum = self.end - length
            values = array("Q")
            if length > 0 and maximum >= self.start:
                for candidate_start in range(self.start, maximum + 1):
                    if self._start_is_valid(candidate_start, length):
                        values.append(candidate_start)
            self._uniform_cache[length] = values
            cached = values
        return cached

    def anchor_candidates(
        self, dinuc: str, length: int, anchor: str, original: Fragment
    ) -> list[int]:
        cached = self._anchor_starts(dinuc, length, anchor)
        output: list[int] = []
        for candidate_start in cached:
            candidate = (int(candidate_start), int(candidate_start) + length)
            if candidate != original and self._coordinate_available(candidate):
                output.append(int(candidate_start))
        return output

    def _anchor_starts(self, dinuc: str, length: int, anchor: str) -> array:
        key = (dinuc, length, anchor)
        cached = self._anchor_cache.get(key)
        if cached is None:
            values = array("Q")
            for relative_position in self._dinuc_index.get(dinuc, ()):
                if anchor == "start":
                    candidate_start = self.start + relative_position
                elif anchor == "end":
                    candidate_start = self.start + relative_position + 2 - length
                else:
                    raise ValueError("anchor must be start or end")
                if self._start_is_valid(candidate_start, length):
                    values.append(candidate_start)
            self._anchor_cache[key] = values
            cached = values
        return cached

    def _choose_start(
        self,
        cached: Sequence[int],
        length: int,
        original: Fragment,
        rng: random.Random,
    ) -> int | None:
        """Choose without copying the cached candidate array per fragment."""
        if not cached:
            return None
        original_start = int(original[0])
        index = bisect_left(cached, original_start)
        original_present = (
            index < len(cached)
            and int(cached[index]) == original_start
            and original[1] - original[0] == length
        )
        if self.max_per_coordinate == 0:
            allowed = len(cached) - int(original_present)
            if allowed <= 0:
                return None
            draw = rng.randrange(allowed)
            if original_present and draw >= index:
                draw += 1
            return int(cached[draw])

        available = [
            int(start)
            for start in cached
            if (int(start), int(start) + length) != original
            and self._coordinate_available((int(start), int(start) + length))
        ]
        return rng.choice(available) if available else None

    def choose_uniform_start(
        self, length: int, original: Fragment, rng: random.Random
    ) -> int | None:
        return self._choose_start(self._uniform_starts(length), length, original, rng)

    def choose_anchor_start(
        self,
        dinuc: str,
        length: int,
        anchor: str,
        original: Fragment,
        rng: random.Random,
    ) -> int | None:
        return self._choose_start(
            self._anchor_starts(dinuc, length, anchor), length, original, rng
        )

    def has_anchor_candidate(
        self, dinuc: str, length: int, anchor: str, original: Fragment
    ) -> bool:
        cached = self._anchor_starts(dinuc, length, anchor)
        if self.max_per_coordinate == 0:
            if not cached:
                return False
            if len(cached) > 1:
                return True
            start = int(cached[0])
            return (start, start + length) != original
        return any(
            (int(start), int(start) + length) != original
            and self._coordinate_available((int(start), int(start) + length))
            for start in cached
        )

    def register(self, fragment: Fragment) -> None:
        assert self.coordinate_counts is not None
        self.coordinate_counts[fragment] = self.coordinate_counts.get(fragment, 0) + 1


def uniform_randomize_fragment(
    fragment: Fragment,
    block: RandomizationBlock,
    rng: random.Random,
) -> Fragment | None:
    length = fragment[1] - fragment[0]
    start = block.choose_uniform_start(length, fragment, rng)
    if start is None:
        return None
    output = (start, start + length)
    block.register(output)
    return output


def place_dinucleotide_matched(
    fragment: Fragment,
    *,
    start_dinuc: str,
    end_dinuc: str,
    block: RandomizationBlock,
    rng: random.Random,
    anchor_prob_start: float = 0.5,
    fallback: str = "uniform",
) -> tuple[Fragment | None, str, str | None, str | None, str | None]:
    """Place one fragment and report status and anchor-selection details.

    Returns ``(fragment, status, selected_anchor, matched_anchor, reason)``.
    ``status`` is ``matched``, ``fallback`` or ``skipped``.
    """
    if fallback not in {"uniform", "skip"}:
        raise ValueError("fallback must be uniform or skip")
    if not 0.0 <= anchor_prob_start <= 1.0:
        raise ValueError("anchor_prob_start must be between 0 and 1")
    length = fragment[1] - fragment[0]
    if length <= 0 or length > block.length:
        return None, "skipped", None, None, "invalid_length"

    start_dinuc = start_dinuc.upper()
    end_dinuc = end_dinuc.upper()
    start_available = block.has_anchor_candidate(start_dinuc, length, "start", fragment)
    end_available = block.has_anchor_candidate(end_dinuc, length, "end", fragment)

    selected = "start" if rng.random() < anchor_prob_start else "end"
    if selected == "start" and start_available:
        chosen_anchor = "start"
    elif selected == "end" and end_available:
        chosen_anchor = "end"
    elif start_available:
        chosen_anchor = "start"
    elif end_available:
        chosen_anchor = "end"
    else:
        if start_dinuc not in DINUCS and end_dinuc not in DINUCS:
            reason = "ambiguous_boundaries"
        else:
            reason = "no_valid_dinucleotide_candidate"
        if fallback == "uniform":
            output = uniform_randomize_fragment(fragment, block, rng)
            return output, "fallback" if output is not None else "skipped", selected, None, reason
        return None, "skipped", selected, None, reason

    dinuc = start_dinuc if chosen_anchor == "start" else end_dinuc
    start = block.choose_anchor_start(
        dinuc, length, chosen_anchor, fragment, rng
    )
    if start is None:  # defensive: availability was checked immediately above
        return None, "skipped", selected, None, "candidate_became_unavailable"
    output = (start, start + length)
    block.register(output)
    return output, "matched", selected, chosen_anchor, None


def uniform_randomize_fragments(
    fragments: Sequence[Fragment], start: int, end: int, *, sequence: str | None = None,
    rng: random.Random | None = None,
) -> List[Fragment]:
    """Compatibility batch helper using canonical sequence when supplied."""
    rng = rng or random
    sequence = sequence if sequence is not None else "A" * (end - start)
    block = RandomizationBlock("region", start, end, sequence)
    output: list[Fragment] = []
    for fragment in fragments:
        randomized = uniform_randomize_fragment(fragment, block, rng)  # type: ignore[arg-type]
        if randomized is not None:
            output.append(randomized)
    return output


def dinuc_anchor_randomize_fragments(
    fragments: Sequence[Fragment],
    start: int,
    end: int,
    window_sequence: str,
    dinuc_positions: Mapping[str, Sequence[int]] | None = None,
    anchor_prob_start: float = 0.5,
    max_anchor_tries: int = 30,
    fallback: str = "uniform",
    *,
    rng: random.Random | None = None,
) -> List[Fragment]:
    """Batch helper using direct valid-candidate selection.

    Candidate positions are filtered before sampling. ``dinuc_positions`` and
    ``max_anchor_tries`` are accepted by this helper but do not alter sampling.
    """
    del dinuc_positions, max_anchor_tries
    rng = rng or random
    block = RandomizationBlock("region", start, end, window_sequence)
    output: list[Fragment] = []
    for fragment in fragments:
        relative_start = fragment[0] - start
        relative_end = fragment[1] - start
        if not (0 <= relative_start <= len(window_sequence) - 2 and 2 <= relative_end <= len(window_sequence)):
            randomized = uniform_randomize_fragment(fragment, block, rng) if fallback == "uniform" else None  # type: ignore[arg-type]
        else:
            randomized, _status, _selected, _matched, _reason = place_dinucleotide_matched(
                fragment,
                start_dinuc=window_sequence[relative_start : relative_start + 2],
                end_dinuc=window_sequence[relative_end - 2 : relative_end],
                block=block,
                rng=rng,  # type: ignore[arg-type]
                anchor_prob_start=anchor_prob_start,
                fallback=fallback,
            )
        if randomized is not None:
            output.append(randomized)
    return output
