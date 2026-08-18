import numpy as np
import pytest

from nucleosuite.percentiles import (
    equal_rank_bins,
    randomized_score_order,
    rank_bins_from_boundaries,
)


def test_equal_rank_bins_cover_all_ranks_with_nearly_equal_counts():
    bins = equal_rank_bins(103, 5)

    assert len(bins) == 20
    assert bins[0].percentile_lower == 0
    assert bins[-1].percentile_upper == 100
    assert bins[0].rank_start == 0
    assert bins[-1].rank_stop == 103
    sizes = [item.rank_stop - item.rank_start for item in bins]
    assert min(sizes) == 5
    assert max(sizes) == 6


def test_equal_rank_bin_size_must_divide_100():
    with pytest.raises(ValueError, match="divide 100"):
        equal_rank_bins(100, 3)


def test_randomized_score_order_is_reproducible_and_randomizes_ties():
    scores = np.ones(30, dtype=float)

    first = randomized_score_order(scores, 11)
    repeated = randomized_score_order(scores, 11)
    changed = randomized_score_order(scores, 12)

    assert np.array_equal(first, repeated)
    assert not np.array_equal(first, np.arange(scores.size))
    assert not np.array_equal(first, changed)


def test_randomized_score_order_preserves_ascending_score_order():
    scores = np.asarray([2, 1, 2, 3, 1, 3, 2], dtype=float)
    order = randomized_score_order(scores, 4)

    assert np.array_equal(scores[order], np.sort(scores))


def test_rank_bins_from_arbitrary_boundaries_use_requested_global_percentages():
    bins = rank_bins_from_boundaries(20, (0, 10, 30, 100))
    assert [(item.percentile_lower, item.percentile_upper) for item in bins] == [
        (0, 10),
        (10, 30),
        (30, 100),
    ]
    assert [(item.rank_start, item.rank_stop) for item in bins] == [
        (0, 2),
        (2, 6),
        (6, 20),
    ]


def test_rank_bins_from_boundaries_allocate_rounding_without_overlap():
    bins = rank_bins_from_boundaries(10, (0, 25, 50, 75, 100))
    assert [item.rank_stop - item.rank_start for item in bins] == [3, 3, 2, 2]
    assert bins[0].rank_start == 0
    assert bins[-1].rank_stop == 10
