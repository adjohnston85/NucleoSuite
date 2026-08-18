from __future__ import annotations

import numpy as np


def test_blacklist_masks_signal_bases_as_missing(tmp_path):
    from nucleosuite.core.blacklist import load_blacklist_unbounded

    path = tmp_path / "blacklist.bed"
    path.write_text("chr1\t2\t5\n")
    index = load_blacklist_unbounded(path)
    values = np.arange(8, dtype=float)
    assert index is not None
    assert index.mask_values("1", 0, values) == 3
    assert np.array_equal(np.flatnonzero(np.isnan(values)), np.array([2, 3, 4]))
    assert np.array_equal(index.valid_mask("chr1", 0, 8), np.array([1, 1, 0, 0, 0, 1, 1, 1], dtype=bool))


def test_dac_masked_opportunities_match_brute_force():
    from nucleosuite.dac import opportunity_vector_from_mask

    valid = np.array([True, False, True, True, False, True])
    got = opportunity_vector_from_mask(valid, 5)
    expected = np.zeros(6)
    for lag in range(1, 6):
        expected[lag] = sum(valid[i] and valid[i + lag] for i in range(len(valid) - lag))
    assert np.array_equal(got, expected)


def test_dcc_masked_opportunities_match_brute_force():
    from nucleosuite.dcc import signed_opportunity_vector_from_masks

    valid_a = np.array([True, False, True, True, False])
    valid_b = np.array([False, True, True, False, True])
    dmax = 4
    got = signed_opportunity_vector_from_masks(valid_a, valid_b, dmax)
    expected = []
    for lag in range(-dmax, dmax + 1):
        expected.append(
            sum(
                bool(valid_a[index] and valid_b[index + lag])
                for index in range(len(valid_a))
                if 0 <= index + lag < len(valid_b)
            )
        )
    assert np.array_equal(got, np.asarray(expected))


def test_bigwig_readers_treat_ordinary_missing_values_as_zero():
    from nucleosuite import dac, dcc

    class Handle:
        def values(self, _chrom, _start, _end, numpy=False):
            values = np.array([1.0, np.nan, np.inf, -np.inf])
            return values if numpy else values.tolist()

    handle = Handle()
    dac_values, _ = dac.read_bigwig_region(
        dac.Track("signal.bw", handle, {"chr1": 4}, {}),
        dac.Region("chr1", 0, 4, "all"),
        None,
    )
    dcc_values, _ = dcc.read_bigwig_region(
        dcc.BigWigTrack("signal.bw", handle, {"chr1": 4}, {}),
        dcc.Region("chr1", 0, 4, "all"),
        None,
    )
    expected = np.array([1.0, 0.0, 0.0, 0.0])
    assert np.array_equal(dac_values, expected)
    assert np.array_equal(dcc_values, expected)


def test_hg19_detection_requires_exact_known_lengths():
    from nucleosuite.core.blacklist import is_hg19_reference

    assert is_hg19_reference({"chr1": 249250621, "chr17": 81195210})
    assert is_hg19_reference({"1": 249250621})
    assert not is_hg19_reference({"chr1": 248956422})
    assert not is_hg19_reference({"chrUn_test": 1000})
