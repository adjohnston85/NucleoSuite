"""Coverage, dyad and fragment-end track calculations."""

from __future__ import annotations

from typing import Dict

import numpy as np

from nucleosuite.sequence.dinucleotide import fragment_dyad

BASIC_TRACKS = (
    "coverage",
    "dyad",
    "fragment_ends",
    "fragment_left_ends",
    "fragment_right_ends",
)
SPARSE_TRACKS = frozenset(
    {"dyad", "fragment_ends", "fragment_left_ends", "fragment_right_ends"}
)


def new_arrays(reference_length: int) -> Dict[str, np.ndarray]:
    return {
        "coverage": np.zeros(reference_length, dtype=np.int32),
        "dyad": np.zeros(reference_length, dtype=np.float64),
        "fragment_ends": np.zeros(reference_length, dtype=np.int32),
        "fragment_left_ends": np.zeros(reference_length, dtype=np.int32),
        "fragment_right_ends": np.zeros(reference_length, dtype=np.int32),
    }


def add_fragment(
    arrays: Dict[str, np.ndarray],
    fragment_start: int,
    fragment_end: int,
    window_start: int,
    window_end: int,
    even_dyad: str = "split",
) -> None:
    """Add a complete fragment to basic tracks.

    Tracks are updated only when the complete fragment lies within the
    adjusted scoring window.
    """
    if even_dyad not in {"split", "left", "right"}:
        raise ValueError("even_dyad must be split, left or right")
    if fragment_start < window_start or fragment_end > window_end:
        return
    length = fragment_end - fragment_start
    if length <= 0:
        return

    left = fragment_start - window_start
    right = fragment_end - window_start
    arrays["coverage"][left:right] += 1

    if length % 2:
        centre = fragment_dyad(fragment_start, fragment_end) - window_start
        if 0 <= centre < len(arrays["dyad"]):
            arrays["dyad"][centre] += 1.0
    else:
        right_centre = fragment_dyad(fragment_start, fragment_end) - window_start
        left_centre = right_centre - 1
        if even_dyad == "split":
            if 0 <= left_centre < len(arrays["dyad"]):
                arrays["dyad"][left_centre] += 0.5
            if 0 <= right_centre < len(arrays["dyad"]):
                arrays["dyad"][right_centre] += 0.5
        else:
            centre = left_centre if even_dyad == "left" else right_centre
            if 0 <= centre < len(arrays["dyad"]):
                arrays["dyad"][centre] += 1.0

    left_end = fragment_start - window_start
    right_end = fragment_end - 1 - window_start
    if 0 <= left_end < len(arrays["fragment_ends"]):
        arrays["fragment_ends"][left_end] += 1
        arrays["fragment_left_ends"][left_end] += 1
    if 0 <= right_end < len(arrays["fragment_ends"]):
        arrays["fragment_ends"][right_end] += 1
        arrays["fragment_right_ends"][right_end] += 1


def cap_sparse_arrays(arrays: Dict[str, np.ndarray], maximum: int) -> None:
    """Cap sparse dyad/end signals in place; ``0`` leaves them unlimited."""
    if maximum < 0:
        raise ValueError("maximum must be 0 or greater")
    if maximum == 0:
        return
    for track in SPARSE_TRACKS:
        np.minimum(arrays[track], maximum, out=arrays[track])


def to_scores(arrays: Dict[str, np.ndarray], contig: str, start: int):
    return {name: [(contig, start, values)] for name, values in arrays.items()}
