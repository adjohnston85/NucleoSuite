import numpy as np

from nucleosuite.scoring.basic_tracks import add_fragment, new_arrays
from nucleosuite.scoring.pns import (
    balanced_boxcar_probability_kernel,
    endpoint_probability_triangle,
    precompute_distributions,
    triangular_probability_kernel,
)
from nucleosuite.scoring.wps import wps_kernel_kircher_exact


def test_even_dyad_split():
    arrays = new_arrays(20)
    add_fragment(arrays, 4, 10, 0, 20, even_dyad="split")
    assert arrays["dyad"][6] == 0.5
    assert arrays["dyad"][7] == 0.5


def test_pns_endpoint_triangle_is_probability_normalized_with_single_odd_midpoint():
    triangle = endpoint_probability_triangle(167)
    assert len(triangle) == 167
    assert np.isclose(np.sum(triangle), 0.5)
    assert np.argmax(triangle) == 83
    assert np.count_nonzero(np.isclose(triangle, triangle.max())) == 1
    assert np.isclose(triangle[82], triangle[84])
    assert triangle[0] == 0.0
    assert triangle[-1] == 0.0


def test_pns_endpoint_triangle_has_two_equal_centres_for_even_mode():
    triangle = endpoint_probability_triangle(168)
    assert np.isclose(np.sum(triangle), 0.5)
    assert np.isclose(triangle[83], triangle[84])
    assert np.count_nonzero(np.isclose(triangle, triangle.max())) == 2


def test_pns_uncentred_kernel_has_unit_probability_mass():
    centred, positive = precompute_distributions([137, 167, 197], 167)
    for length in (137, 167, 197):
        assert np.isclose(np.sum(positive[length]), 1.0)
        assert np.isclose(np.sum(centred[length]), 0.0)


def test_pns_kernels_are_mean_centred():
    centred, positive = precompute_distributions([167], 167)
    assert np.isclose(np.mean(centred[167]), 0.0)
    assert len(positive[167]) == 167


def test_wps_kernel_geometry():
    kernel = wps_kernel_kircher_exact(167, protection=120)
    assert len(kernel) == 286
    assert (kernel == 1).sum() == 48
    assert (kernel == -1).sum() == 238


def test_bns_uncentred_kernel_has_unit_mass_and_centred_kernel_sums_to_zero():
    centred, positive = precompute_distributions([120, 167, 180], 167, "bns")
    assert len(positive[120]) == 214
    assert len(positive[167]) == 167
    assert len(positive[180]) == 180
    for length in (120, 167, 180):
        assert np.isclose(np.sum(positive[length]), 1.0)
        assert np.isclose(np.sum(centred[length]), 0.0)
        assert np.allclose(positive[length], positive[length][::-1])
        assert np.allclose(centred[length], centred[length][::-1])


def test_bns_120_bp_mode167_has_symmetric_balanced_geometry():
    centred, positive = precompute_distributions([120], 167, "bns")
    kernel = centred[120]
    assert len(kernel) == 214
    assert np.sum(kernel > 0) == 106
    assert np.sum(kernel < 0) == 106
    assert np.sum(np.isclose(kernel, 0.0)) == 2
    assert np.all(kernel[:53] < 0)
    assert np.isclose(kernel[53], 0.0)
    assert np.all(kernel[54:160] > 0)
    assert np.isclose(kernel[160], 0.0)
    assert np.all(kernel[161:] < 0)
    assert np.isclose(np.sum(positive[120]), 1.0)


def test_bns_even_support_divisible_by_four_has_equal_positive_and_negative_counts():
    raw = balanced_boxcar_probability_kernel(180)
    centred = raw - raw.mean()
    assert np.sum(centred > 0) == 90
    assert np.sum(centred < 0) == 90
    assert np.all(centred[:45] < 0)
    assert np.all(centred[45:135] > 0)
    assert np.all(centred[135:] < 0)


def test_bns_odd_support_uses_symmetric_half_weight_boundary_values():
    expected = {
        167: (-0.5, 83, 84),
        169: (0.5, 85, 84),
        177: (0.5, 89, 88),
    }
    for support, (edge_multiple, positive_count, negative_count) in expected.items():
        raw = balanced_boxcar_probability_kernel(support)
        centred = raw - raw.mean()
        amplitude = 1.0 / support
        assert np.allclose(centred, centred[::-1])
        assert np.isclose(np.sum(centred), 0.0)
        assert np.sum(centred > 0) == positive_count
        assert np.sum(centred < 0) == negative_count
        partial = np.flatnonzero(
            np.isclose(np.abs(centred), 0.5 * amplitude)
        )
        assert len(partial) == 2
        assert np.allclose(centred[partial], edge_multiple * amplitude)


def test_bns_178_bp_support_has_two_zero_boundaries_and_balanced_full_values():
    raw = balanced_boxcar_probability_kernel(178)
    centred = raw - raw.mean()
    assert np.sum(centred > 0) == 88
    assert np.sum(centred < 0) == 88
    assert np.sum(np.isclose(centred, 0.0)) == 2
    assert np.isclose(np.sum(centred), 0.0)


def test_tns_triangle_has_unit_mass_and_zero_boundaries():
    for support in (167, 180, 197, 214):
        raw = triangular_probability_kernel(support)
        assert len(raw) == support
        assert np.isclose(np.sum(raw), 1.0)
        assert np.isclose(raw[0], 0.0)
        assert np.isclose(raw[-1], 0.0)
        assert np.allclose(raw, raw[::-1])


def test_tns_odd_support_has_one_central_maximum():
    raw = triangular_probability_kernel(167)
    assert np.argmax(raw) == 83
    assert np.count_nonzero(np.isclose(raw, raw.max())) == 1


def test_tns_even_support_has_two_base_central_plateau():
    raw = triangular_probability_kernel(180)
    assert np.isclose(raw[89], raw[90])
    assert np.count_nonzero(np.isclose(raw, raw.max())) == 2


def test_tns_precomputed_kernels_are_unit_mass_and_mean_centred():
    centred, positive = precompute_distributions([120, 167, 180], 167, "tns")
    assert len(positive[120]) == 214
    assert len(positive[167]) == 167
    assert len(positive[180]) == 180
    for length in (120, 167, 180):
        assert np.isclose(np.sum(positive[length]), 1.0)
        assert np.isclose(np.sum(centred[length]), 0.0)
        assert np.allclose(positive[length], positive[length][::-1])
        assert np.allclose(centred[length], centred[length][::-1])


def test_tns_137_and_197_mode167_use_identical_197_bp_triangle():
    centred, positive = precompute_distributions([137, 197], 167, "tns")
    assert len(positive[137]) == 197
    assert len(positive[197]) == 197
    assert np.allclose(positive[137], positive[197])
    assert np.allclose(centred[137], centred[197])
