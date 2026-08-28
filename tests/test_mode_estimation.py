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


def test_mode_histogram_is_unsmoothed_by_default_and_smoothing_is_explicit():
    counts = np.asarray([0, 10, 0, 9, 9, 9, 0], dtype=int)
    raw_mode, *_ = bootstrap_histogram_mode(
        counts, lower=120, replicates=20, seed=7
    )
    smoothed_mode, *_ = bootstrap_histogram_mode(
        counts,
        lower=120,
        replicates=20,
        seed=7,
        histogram_smoothing="binomial",
    )

    assert raw_mode == 121
    assert smoothed_mode == 124


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


def test_nuc_score_auto_mode_can_use_fragment_interval_input(tmp_path):
    from nucleosuite.cli.main import build_parser
    from nucleosuite.cli.nuc_score import _resolve_mode_and_fragment_range

    fragments = tmp_path / "fragments.bed"
    fragments.write_text(
        "".join(
            f"chr1\t{index * 200}\t{index * 200 + length}\n"
            for index, length in enumerate([167] * 8 + [166] * 2)
        ),
        encoding="utf-8",
    )
    args = build_parser().parse_args(
        [
            "nuc-score",
            "--fragments",
            str(fragments),
            "--mode-min-fragments",
            "1",
            "--mode-batch-fragments",
            "1",
            "--mode-max-fragments",
            "20",
            "--mode-bootstrap",
            "10",
            "--mode-stable-checkpoints",
            "1",
        ]
    )

    mode, estimate, source, _seed = _resolve_mode_and_fragment_range(args)

    assert mode == 167
    assert source == "automatic"
    assert estimate is not None
    assert estimate.histogram_smoothing == "none"
    assert (args.frag_lower, args.frag_upper) == (137, 197)
