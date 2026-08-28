import numpy as np

from nucleosuite.scoring.basic_tracks import add_fragment as add_basic_fragment
from nucleosuite.scoring.basic_tracks import new_arrays as new_basic_arrays
from nucleosuite.scoring.pns import (
    add_fragment as add_pns_fragment,
    new_arrays as new_pns_arrays,
    pns_nucleosome_kernel,
    precompute_distributions,
    scoring_support_length,
)
from nucleosuite.scoring.wps import wps_kernel_kircher_exact


def test_even_dyad_split():
    arrays = new_basic_arrays(20)
    add_basic_fragment(arrays, 4, 10, 0, 20, even_dyad="split")
    assert arrays["dyad"][6] == 0.5
    assert arrays["dyad"][7] == 0.5


def test_pns_kernel_has_percent_positive_mass_and_balanced_negative_mass():
    for support in (167, 180, 214):
        kernel = pns_nucleosome_kernel(support)
        assert len(kernel) == support
        assert np.allclose(kernel, kernel[::-1])
        assert np.isclose(kernel[kernel > 0].sum(), 100.0)
        assert np.isclose(kernel[kernel < 0].sum(), -100.0)
        assert np.isclose(np.abs(kernel).sum(), 200.0)
        assert np.isclose(kernel.sum(), 0.0)
        assert kernel[0] < 0
        assert np.isclose(kernel[0], kernel[-1])


def test_pns_precomputed_distributions_use_the_fixed_kernel():
    centred, positive = precompute_distributions([120, 167, 180], 167)
    assert len(centred[120]) == 214
    assert len(centred[167]) == 167
    assert len(centred[180]) == 180
    for length in (120, 167, 180):
        assert np.isclose(centred[length].sum(), 0.0)
        assert np.isclose(centred[length][centred[length] > 0].sum(), 100.0)
        assert np.isclose(centred[length][centred[length] < 0].sum(), -100.0)
        assert np.allclose(centred[length], centred[length][::-1])
        assert np.allclose(positive[length], positive[length][::-1])


def test_pospns_is_a_nonnegative_vertical_shift_without_renormalization():
    centred, positive = precompute_distributions([167], 167)
    signed = centred[167]
    shifted = positive[167]
    assert np.allclose(shifted - signed, -signed.min())
    assert np.isclose(shifted.min(), 0.0)
    assert np.all(shifted >= 0)
    assert not np.isclose(shifted.sum(), 100.0)


def test_pns_support_broadens_for_fragments_away_from_mode():
    assert scoring_support_length(120, 167) == 214
    assert scoring_support_length(167, 167) == 167
    assert scoring_support_length(180, 167) == 180
    for length in range(20, 221):
        assert scoring_support_length(length, 167) % 2 == length % 2


def test_pns_short_fragment_extends_symmetrically_beyond_fragment():
    centred, positive = precompute_distributions([120], 167)
    arrays = new_pns_arrays(400)
    add_pns_fragment(
        arrays, 100, 220, 0, 400, 167, centred, positive
    )
    nonzero = np.flatnonzero(arrays["pns"] != 0)
    assert nonzero[0] == 53
    assert nonzero[-1] == 266
    assert np.allclose(arrays["pns"][53:267], arrays["pns"][53:267][::-1])


def test_wps_kernel_geometry():
    kernel = wps_kernel_kircher_exact(167, protection=120)
    assert len(kernel) == 286
    assert (kernel == 1).sum() == 48
    assert (kernel == -1).sum() == 238
