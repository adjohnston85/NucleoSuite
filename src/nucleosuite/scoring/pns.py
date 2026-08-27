"""Nucleosome scoring kernels and array operations for PNS, SNS, BNS, and TNS."""

from __future__ import annotations

from typing import Dict, Iterable, Tuple

import numpy as np
from scipy.signal import savgol_filter

from nucleosuite.sequence.dinucleotide import fragment_dyad

PNS_TRACKS = ("pns_smoothed", "pns", "posPNS")
SNS_TRACKS = ("sns_smoothed", "sns", "posSNS")
BNS_TRACKS = ("bns_smoothed", "bns", "posBNS")
TNS_TRACKS = ("tns_smoothed", "tns", "posTNS")
SCORING_METHODS = ("sns", "pns", "bns", "tns")


def scoring_track_names(scoring_method: str) -> tuple[str, str, str]:
    """Return ``(smoothed, score, positive-reference)`` track names for a method."""
    if scoring_method == "pns":
        return PNS_TRACKS
    if scoring_method == "sns":
        return SNS_TRACKS
    if scoring_method == "bns":
        return BNS_TRACKS
    if scoring_method == "tns":
        return TNS_TRACKS
    raise ValueError(f"Unknown nucleosome scoring method: {scoring_method}")


def scoring_support_length(fragment_length: int, mode_dna_length: int) -> int:
    """Return the genomic support used by PNS/SNS/BNS/TNS for one fragment length."""
    return (
        2 * mode_dna_length - fragment_length
        if fragment_length < mode_dna_length
        else fragment_length
    )


def endpoint_probability_triangle(
    mode_dna_length: int, total_length: int | None = None
) -> np.ndarray:
    """Return one endpoint-derived PNS probability triangle.

    The non-zero triangle spans ``mode_dna_length`` bases.  Its discrete
    probability mass is exactly 0.5, so the two mirrored endpoint triangles
    used for one fragment sum to a total mass of 1 before mean centring.

    Odd mode lengths have one unique central maximum.  Even mode lengths have
    two equal central maxima, because their geometric centre lies between two
    bases.
    """
    if mode_dna_length < 3:
        raise ValueError("mode_dna_length must be at least 3")
    if total_length is None:
        total_length = mode_dna_length
    if total_length < mode_dna_length:
        raise ValueError("total_length must be at least mode_dna_length")

    midpoint = (mode_dna_length - 1) // 2
    half_mass_count = mode_dna_length // 2
    probability_denominator = 2.0 * midpoint * half_mass_count

    triangle = np.zeros(total_length, dtype=np.float64)
    for index in range(mode_dna_length):
        triangle[index] = (
            min(index, mode_dna_length - 1 - index) / probability_denominator
        )
    return triangle


def balanced_boxcar_probability_kernel(total_length: int) -> np.ndarray:
    """Return the symmetric unit-mass uncentred BNS boxcar kernel.

    BNS divides the scoring support into equal effective central and flank
    widths.  The raw central kernel has total mass 1 and the outer flanks are
    zero.  Discrete boundary weights keep the kernel centred when the support
    cannot be divided into four equal integer blocks.  Subtracting the kernel
    mean therefore produces balanced positive and negative contributions.
    """
    if total_length < 2:
        raise ValueError("total_length must be at least 2")

    n = int(total_length)
    centred = np.empty(n, dtype=np.float64)
    amplitude = 1.0 / n
    remainder = n % 4
    quarter = n // 4

    if remainder == 0:
        # k negative | 2k positive | k negative
        centred[:quarter] = -amplitude
        centred[quarter : 3 * quarter] = amplitude
        centred[3 * quarter :] = -amplitude
    elif remainder == 1:
        # k negative | half-positive | (2k-1) positive | half-positive | k negative
        centred[:quarter] = -amplitude
        left_edge = quarter
        right_edge = n - quarter - 1
        centred[left_edge] = 0.5 * amplitude
        centred[left_edge + 1 : right_edge] = amplitude
        centred[right_edge] = 0.5 * amplitude
        centred[right_edge + 1 :] = -amplitude
    elif remainder == 2:
        # k negative | zero | 2k positive | zero | k negative
        centred[:quarter] = -amplitude
        left_edge = quarter
        right_edge = n - quarter - 1
        centred[left_edge] = 0.0
        centred[left_edge + 1 : right_edge] = amplitude
        centred[right_edge] = 0.0
        centred[right_edge + 1 :] = -amplitude
    else:
        # k negative | half-negative | (2k+1) positive | half-negative | k negative
        centred[:quarter] = -amplitude
        left_edge = quarter
        right_edge = n - quarter - 1
        centred[left_edge] = -0.5 * amplitude
        centred[left_edge + 1 : right_edge] = amplitude
        centred[right_edge] = -0.5 * amplitude
        centred[right_edge + 1 :] = -amplitude

    # Adding 1/n converts the balanced, zero-sum contribution into a
    # non-negative unit-mass boxcar with zero outer flanks.
    uncentred = centred + amplitude
    if not np.isclose(np.sum(uncentred), 1.0):  # pragma: no cover - invariant guard
        raise RuntimeError("BNS unit-mass normalization failed")
    return uncentred


def triangular_probability_kernel(total_length: int) -> np.ndarray:
    """Return a symmetric unit-mass TNS triangle over ``total_length`` bases.

    The triangle is zero at both support boundaries.  Odd support lengths have
    one central maximum; even support lengths have a two-base central plateau.
    The discrete values are normalized so that the complete uncentred kernel
    has total mass 1 before mean centring.
    """
    if total_length < 3:
        raise ValueError("total_length must be at least 3")

    n = int(total_length)
    positions = np.arange(n, dtype=np.float64)
    triangle = np.minimum(positions, positions[::-1])
    total_mass = float(np.sum(triangle))
    if total_mass <= 0:  # pragma: no cover - guarded by total_length >= 3
        raise RuntimeError("TNS triangle normalization failed")
    triangle /= total_mass
    return triangle


def sinusoidal_nucleosome_kernel(total_length: int) -> np.ndarray:
    """Return one discrete, symmetric, zero-sum SNS kernel.

    The raw shape is one inverted cosine cycle sampled at integer genomic bins,
    with its minimum at both support boundaries and its maximum at the centre.
    Positive and negative samples are normalized separately so that the
    positive mass is exactly +50 and the negative mass is exactly -50.
    Consequently every complete fragment contributes zero net mass and total
    absolute mass 100, independent of support width.
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
        raise RuntimeError("SNS signed-mass normalization failed")

    kernel = np.zeros(n, dtype=np.float64)
    kernel[positive_mask] = raw[positive_mask] * (50.0 / positive_mass)
    kernel[negative_mask] = raw[negative_mask] * (50.0 / negative_mass)

    if not np.allclose(kernel, kernel[::-1], atol=1e-15, rtol=1e-12):  # pragma: no cover
        raise RuntimeError("SNS symmetry normalization failed")
    if not np.isclose(np.sum(kernel), 0.0, atol=1e-12):  # pragma: no cover
        raise RuntimeError("SNS zero-sum normalization failed")
    return kernel


def precompute_distributions(
    fragment_lengths: Iterable[int],
    mode_dna_length: int,
    scoring_method: str = "sns",
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
        if scoring_method == "pns":
            left_triangle = endpoint_probability_triangle(mode_dna_length, total_length)
            combined = left_triangle + left_triangle[::-1]
            positive[fragment_length] = combined.copy()
            centred[fragment_length] = combined - np.mean(combined)
        elif scoring_method == "sns":
            signed = sinusoidal_nucleosome_kernel(total_length)
            # posSNS preserves the complete SNS waveform and shifts it upward
            # by the per-fragment minimum so that every value is non-negative.
            # It is deliberately not clipped and is not renormalized after the
            # vertical translation.
            positive_reference = signed - float(np.min(signed))
            positive[fragment_length] = positive_reference
            centred[fragment_length] = signed
        elif scoring_method == "bns":
            combined = balanced_boxcar_probability_kernel(total_length)
            positive[fragment_length] = combined.copy()
            centred[fragment_length] = combined - np.mean(combined)
        else:
            combined = triangular_probability_kernel(total_length)
            positive[fragment_length] = combined.copy()
            centred[fragment_length] = combined - np.mean(combined)
    return centred, positive


def new_arrays(reference_length: int, scoring_method: str = "sns"):
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
    scoring_method: str = "sns",
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
    scoring_method: str = "sns",
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
