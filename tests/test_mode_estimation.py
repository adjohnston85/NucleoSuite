import numpy as np

from nucleosuite.mode_estimation import (
    bootstrap_histogram_mode,
    estimate_mode_from_lengths,
    pooled_mode_estimate,
)


def test_bootstrap_mode_recovers_dominant_length():
    counts = np.zeros(131, dtype=int)
    counts[47] = 1000  # 120 + 47 = 167
    counts[46] = 500
    counts[48] = 450
    mode, low, high, modes = bootstrap_histogram_mode(
        counts, lower=120, replicates=100, seed=7
    )
    assert mode == 167
    assert low <= 167 <= high
    assert modes.size == 100


def test_pooled_mode_equal_weights_target_and_control_histograms():
    target = estimate_mode_from_lengths(
        [166] * 1000 + [167] * 500,
        bootstrap_replicates=50,
        seed=3,
    )
    control = estimate_mode_from_lengths(
        [168] * 100 + [167] * 50,
        bootstrap_replicates=50,
        seed=4,
    )
    pooled = pooled_mode_estimate(target, control, bootstrap_replicates=50, seed=5)

    assert target.mode == 166
    assert control.mode == 168
    assert 166 <= pooled.mode <= 168
    # Equal histogram weighting prevents the deeper target from dominating by depth.
    assert pooled.mode != 166 or pooled.mode_search_fragments <= target.mode_search_fragments
