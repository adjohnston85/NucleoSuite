"""Length-adaptive probabilistic nucleosome scoring."""

from __future__ import annotations

from typing import Dict, Iterable, Tuple

import numpy as np
from scipy.signal import savgol_filter

from nucleosuite.sequence.dinucleotide import fragment_dyad

PNS_TRACKS = ("pns_smoothed", "pns", "posPNS")
SCORING_METHODS = ("pns",)


def scoring_track_names(scoring_method: str) -> tuple[str, str, str]:
    """Return ``(smoothed, score, positive-reference)`` PNS track names."""
    if scoring_method != "pns":
        raise ValueError(f"Unknown nucleosome scoring method: {scoring_method}")
    return PNS_TRACKS


def scoring_support_length(fragment_length: int, mode_dna_length: int) -> int:
    """Return the genomic support used by PNS for one fragment length."""
    return (
        2 * mode_dna_length - fragment_length
        if fragment_length < mode_dna_length
        else fragment_length
    )


def pns_nucleosome_kernel(total_length: int) -> np.ndarray:
    """Return one discrete, symmetric, zero-sum PNS kernel.

    The raw shape is one inverted cosine cycle sampled at integer genomic bins,
    with its minimum at both support boundaries and its maximum at the centre.
    Positive and negative samples are normalized separately so that the
    positive mass is exactly +100 and the negative mass is exactly -100.
    The positive distribution therefore represents probability in percent.
    Every complete fragment contributes zero net mass and total absolute mass
    200, independent of support width.
    """
    if total_length < 3:
        raise ValueError("total_length must be at least 3")

    n = int(total_length)
    positions = np.arange(n, dtype=np.float64)
    raw = -np.cos((2.0 * np.pi * positions) / float(n - 1))
    raw[np.isclose(raw, 0.0, atol=1e-15, rtol=0.0)] = 0.0

    positive_mask = raw > 0.0
    negative_mask = raw < 0.0
    positive_mass = float(np.sum(raw[positive_mask], dtype=np.float64))
    negative_mass = float(np.sum(-raw[negative_mask], dtype=np.float64))
    if positive_mass <= 0.0 or negative_mass <= 0.0:  # pragma: no cover
        raise RuntimeError("PNS signed-mass normalization failed")

    kernel = np.zeros(n, dtype=np.float64)
    kernel[positive_mask] = raw[positive_mask] * (100.0 / positive_mass)
    kernel[negative_mask] = raw[negative_mask] * (100.0 / negative_mass)

    if not np.allclose(kernel, kernel[::-1], atol=1e-15, rtol=1e-12):  # pragma: no cover
        raise RuntimeError("PNS symmetry normalization failed")
    if not np.isclose(np.sum(kernel), 0.0, atol=1e-12):  # pragma: no cover
        raise RuntimeError("PNS zero-sum normalization failed")
    return kernel


def precompute_distributions(
    fragment_lengths: Iterable[int],
    mode_dna_length: int,
    scoring_method: str = "pns",
) -> Tuple[Dict[int, np.ndarray], Dict[int, np.ndarray]]:
    """Return score and non-negative reference kernels by fragment length."""
    if mode_dna_length < 3:
        raise ValueError("mode_dna_length must be at least 3")
    if scoring_method not in SCORING_METHODS:
        raise ValueError(f"Unknown nucleosome scoring method: {scoring_method}")

    centred: Dict[int, np.ndarray] = {}
    positive: Dict[int, np.ndarray] = {}

    for fragment_length in fragment_lengths:
        total_length = scoring_support_length(fragment_length, mode_dna_length)
        signed = pns_nucleosome_kernel(total_length)
        # posPNS preserves the complete waveform but shifts it upward so that
        # every value is non-negative. It is an auxiliary probability-like
        # reference track and is intentionally not renormalized after shifting.
        positive[fragment_length] = signed - float(np.min(signed))
        centred[fragment_length] = signed
    return centred, positive


def new_arrays(reference_length: int, scoring_method: str = "pns"):
    _, score_track, positive_track = scoring_track_names(scoring_method)
    return {
        score_track: np.zeros(reference_length, dtype=np.float64),
        positive_track: np.zeros(reference_length, dtype=np.float64),
    }


def add_fragment(
    arrays,
    fragment_start: int,
    fragment_end: int,
    window_start: int,
    window_end: int,
    mode_dna_length: int,
    centred_distributions,
    positive_distributions,
    scoring_method: str = "pns",
) -> None:
    fragment_length = fragment_end - fragment_start
    centred_kernel = centred_distributions.get(fragment_length)
    positive_kernel = positive_distributions.get(fragment_length)
    if centred_kernel is None or positive_kernel is None:
        return

    _, score_track, positive_track = scoring_track_names(scoring_method)
    reference_length = window_end - window_start
    if fragment_length < mode_dna_length:
        total_length = scoring_support_length(fragment_length, mode_dna_length)
        centre = fragment_dyad(fragment_start, fragment_end) - window_start
        start_position = centre - total_length // 2
        end_position = start_position + total_length
    else:
        start_position = fragment_start - window_start
        end_position = fragment_end - window_start

    centred_values = centred_kernel
    positive_values = positive_kernel
    if start_position < 0:
        trim = -start_position
        centred_values = centred_values[trim:]
        positive_values = positive_values[trim:]
        start_position = 0
    if end_position > reference_length:
        trim_length = max(0, reference_length - start_position)
        centred_values = centred_values[:trim_length]
        positive_values = positive_values[:trim_length]

    if 0 <= start_position < reference_length and centred_values.size:
        stop = start_position + centred_values.size
        arrays[score_track][start_position:stop] += centred_values
        arrays[positive_track][start_position:stop] += positive_values


def to_scores(
    arrays,
    contig: str,
    start: int,
    smooth_window: int = 0,
    smooth_order: int = 2,
    scoring_method: str = "pns",
):
    smoothed_track, score_track, positive_track = scoring_track_names(scoring_method)
    score_values = arrays[score_track]
    if (
        smooth_window >= 3
        and smooth_window % 2 == 1
        and len(score_values) >= smooth_window
        and smooth_order < smooth_window
    ):
        smoothed = savgol_filter(score_values, smooth_window, smooth_order)
    else:
        smoothed = score_values.copy()
    return {
        smoothed_track: [(contig, start, smoothed)],
        score_track: [(contig, start, score_values)],
        positive_track: [(contig, start, arrays[positive_track])],
    }
