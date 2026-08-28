import numpy as np

from nucleosuite.scoring.basic_tracks import add_fragment, new_arrays
from nucleosuite.scoring.pns import (
    precompute_distributions,
    sinusoidal_nucleosome_kernel,
)
from nucleosuite.scoring.wps import wps_kernel_kircher_exact


def test_even_dyad_split():
    arrays = new_arrays(20)
    add_fragment(arrays, 4, 10, 0, 20, even_dyad="split")
    assert arrays["dyad"][6] == 0.5
    assert arrays["dyad"][7] == 0.5










def test_wps_kernel_geometry():
    kernel = wps_kernel_kircher_exact(167, protection=120)
    assert len(kernel) == 286
    assert (kernel == 1).sum() == 48
    assert (kernel == -1).sum() == 238






















def test_pns_discrete_kernel_has_exact_balanced_signed_mass():
    for support in range(3, 1001):
        kernel = sinusoidal_nucleosome_kernel(support)
        assert len(kernel) == support
        assert np.allclose(kernel, kernel[::-1])
        assert np.isclose(kernel[kernel > 0].sum(), 100.0)
        assert np.isclose(kernel[kernel < 0].sum(), -100.0)
        assert np.isclose(np.abs(kernel).sum(), 200.0)
        assert np.isclose(kernel.sum(), 0.0)
        assert kernel[0] < 0
        assert np.isclose(kernel[0], kernel[-1])


def test_pns_mode167_support_broadens_on_both_sides_of_167():
    centred, positive = precompute_distributions([120, 167, 180], 167, "pns")
    assert len(centred[120]) == 214
    assert len(centred[167]) == 167
    assert len(centred[180]) == 180
    for length in (120, 167, 180):
        assert np.allclose(centred[length], centred[length][::-1])
        assert np.isclose(centred[length][centred[length] > 0].sum(), 100.0)
        assert np.isclose(centred[length][centred[length] < 0].sum(), -100.0)
        assert np.allclose(positive[length], centred[length] - centred[length].min())
        assert np.isclose(positive[length].min(), 0.0)
        assert np.all(positive[length] >= 0)
        assert np.allclose(positive[length], positive[length][::-1])


def test_pospns_is_vertical_shift_of_complete_signed_wave_without_clipping_or_renormalization():
    centred, positive = precompute_distributions([167], 167, "pns")
    signed = centred[167]
    shifted = positive[167]
    assert np.allclose(shifted - signed, -signed.min())
    assert np.isclose(shifted.min(), 0.0)
    assert np.count_nonzero(shifted > 0) > np.count_nonzero(signed > 0)
    assert not np.isclose(shifted.sum(), 1.0)


def test_pns_odd_support_has_one_central_maximum_and_even_support_two():
    odd = sinusoidal_nucleosome_kernel(167)
    even = sinusoidal_nucleosome_kernel(180)
    assert np.argmax(odd) == 83
    assert np.count_nonzero(np.isclose(odd, odd.max())) == 1
    assert np.isclose(even[89], even[90])
    assert np.count_nonzero(np.isclose(even, even.max())) == 2


def test_pns_support_parity_matches_fragment_length_for_exact_midpoint_symmetry():
    from nucleosuite.scoring.pns import scoring_support_length

    for length in range(20, 221):
        support = scoring_support_length(length, 167)
        assert support % 2 == length % 2


def test_pns_short_fragment_extends_symmetrically_beyond_fragment():
    from nucleosuite.scoring.pns import add_fragment as add_score_fragment
    from nucleosuite.scoring.pns import new_arrays as new_score_arrays

    centred, positive = precompute_distributions([120], 167, "pns")
    arrays = new_score_arrays(400, "pns")
    add_score_fragment(
        arrays, 100, 220, 0, 400, 167, centred, positive, scoring_method="pns"
    )
    nonzero = np.flatnonzero(arrays["pns"] != 0)
    assert nonzero[0] == 53
    assert nonzero[-1] == 266
    assert np.allclose(arrays["pns"][53:267], arrays["pns"][53:267][::-1])
