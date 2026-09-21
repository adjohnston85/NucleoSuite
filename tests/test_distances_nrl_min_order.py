"""Regression-order filtering for pooled and interval-contained category NRLs."""
from __future__ import annotations

import csv
from collections import Counter

import pytest

from nucleosuite import distances
from test_distances_nrl_regression import synthetic_results
from test_distances_state_isolated_orders import _write_example


def test_pooled_regression_starts_at_order_two_without_renumbering():
    results = synthetic_results()
    # The +1 peak is intentionally inconsistent with the +2/+3 repeat length.
    results.genome_all[1] = Counter({215: 25})
    unfiltered = distances.collect_nrl_regressions(
        results, max_order=3, include_chromosomes=False,
        include_genome=True, nrl_mode='raw',
    )[0]
    filtered = distances.collect_nrl_regressions(
        results, max_order=3, include_chromosomes=False,
        include_genome=True, nrl_mode='raw', nrl_min_order=2,
    )[0]
    assert [peak.order for peak in filtered.peaks] == [2, 3]
    assert [peak.peak_distance for peak in filtered.peaks] == [370, 555]
    assert filtered.slope == pytest.approx(185)
    assert unfiltered.slope != pytest.approx(filtered.slope)
    assert results.genome_all[1][215] == 25


def test_state_regression_skips_first_order_and_keeps_original_numbers():
    results = synthetic_results()
    results.genome_state = {
        1: {'active_genes': Counter({250: 100})},
        2: {'active_genes': Counter({372: 90})},
        3: {'active_genes': Counter({558: 80})},
    }
    rows = distances.collect_state_nrl_regressions(
        results, state='active_genes', max_order=3,
        include_chromosomes=False, include_genome=True, nrl_mode='raw',
        count_smooth_window=21, count_smooth_polyorder=2,
        min_distance=1, max_distance=1500, nrl_min_order=2,
    )
    assert len(rows) == 1
    assert [peak.order for peak in rows[0].peaks] == [2, 3]
    assert rows[0].slope == pytest.approx(186)


def test_cli_keeps_all_distributions_but_regresses_orders_two_to_three(tmp_path):
    peaks, states = _write_example(tmp_path)
    prefix = tmp_path / 'spacing'
    assert distances.main([
        str(peaks), '--state-bed', str(states), '--position-column', '7',
        '--max-order', '3', '--nrl-min-order', '2', '--max-distance', '1000',
        '--count-smooth-window', '0', '--percent-smooth-window', '0',
        '--nrl-mode', 'raw', '--scope', 'combined_chromosomes',
        '--output-prefix', str(prefix),
    ]) == 0
    pooled_summary = next(tmp_path.glob('*_scorepct0_summary.tsv'))
    with pooled_summary.open() as handle:
        summary = list(csv.DictReader(handle, delimiter='\t'))
    assert {int(row['order']) for row in summary if row['state'] == 'All'} == {1, 2, 3}
    root = next(tmp_path.glob('*_scorepct0_states'))
    for category, expected in [('active_genes', 186), ('repressed_genes', 191)]:
        folder = root / category
        with (folder / 'summary.tsv').open() as handle:
            category_distances = list(csv.DictReader(handle, delimiter='\t'))
        assert {int(row['order']) for row in category_distances} == {1, 2, 3}
        with (folder / 'state_nrl_regression_combined_chromosomes.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        assert [int(row['order']) for row in rows] == [2, 3]
        with (folder / 'state_nrl_regression_summary.tsv').open() as handle:
            fits = list(csv.DictReader(handle, delimiter='\t'))
        assert float(fits[0]['nrl_bp']) == pytest.approx(expected)
        assert int(fits[0]['min_order']) == 2
        assert int(fits[0]['max_order']) == 3
    with next(tmp_path.glob('*_scorepct0_metadata.tsv')).open() as handle:
        metadata = dict(csv.reader(handle, delimiter='\t'))
    assert metadata['nrl_min_order'] == '2'


@pytest.mark.parametrize('minimum,maximum', [(0, 7), (8, 7)])
def test_invalid_nrl_min_order_rejected(tmp_path, minimum, maximum):
    peaks, states = _write_example(tmp_path)
    with pytest.raises(ValueError, match='--nrl-min-order'):
        distances.collect_nrl_regressions(
            synthetic_results(), max_order=maximum,
            include_chromosomes=False, include_genome=True,
            nrl_min_order=minimum,
        )
