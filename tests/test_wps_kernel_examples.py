from __future__ import annotations

import numpy as np

from nucleosuite.scoring.wps import wps_kernel_kircher_exact


def test_wps_kernel_example_lengths_match_documented_shapes():
    k120 = wps_kernel_kircher_exact(120, protection=120)
    k167 = wps_kernel_kircher_exact(167, protection=120)
    k180 = wps_kernel_kircher_exact(180, protection=120)

    assert len(k120) == 239
    assert len(k167) == 286
    assert len(k180) == 299

    assert (k120 == 1).sum() == 1
    assert (k167 == 1).sum() == 48
    assert (k180 == 1).sum() == 61

    assert (k120 == -1).sum() == 238
    assert (k167 == -1).sum() == 238
    assert (k180 == -1).sum() == 238


def test_wps_kernel_is_symmetric_about_fragment_midpoint():
    for length in (120, 167, 180):
        kernel = wps_kernel_kircher_exact(length, protection=120)
        np.testing.assert_array_equal(kernel, kernel[::-1])
