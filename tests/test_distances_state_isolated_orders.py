"""Chromatin-state interval isolation, per-category plots and regressions."""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import pytest

from nucleosuite import distances


def _write_example(tmp_path: Path, *, overlapping: bool = False):
    states = tmp_path / 'states.bed'
    rows = [
        'chr1\t0\t800\tactive_genes\t0\t.\t0\t800\t255,0,0',
        'chr1\t800\t1550\tactive_genes\t0\t.\t800\t1550\t255,0,0',
        'chr1\t1550\t2600\trepressed_genes\t0\t.\t1550\t2600\t0,0,255',
    ]
    if overlapping:
        rows.insert(1, 'chr1\t200\t750\tactive_genes\t0\t.\t200\t750\t255,0,0')
    states.write_text('\n'.join(rows) + '\n')
    centres = [100, 286, 472, 658, 900, 1086, 1272, 1458,
               1650, 1841, 2032, 2223, 2414]
    peaks = tmp_path / 'peaks.bed'
    peaks.write_text(''.join(
        f'1\t{x}\t{x+1}\tp{i}\t10\t.\t{x}\t{x+1}\n'
        for i, x in enumerate(centres)
    ))
    return peaks, states


def test_state_isolated_all_orders_and_independent_nrls(tmp_path):
    peaks, states = _write_example(tmp_path)
    prefix = tmp_path / 'spacing'
    assert distances.main([
        str(peaks), '--state-bed', str(states), '--position-column', '7',
        '--max-order', '3', '--max-distance', '1000',
        '--count-smooth-window', '0', '--percent-smooth-window', '0',
        '--nrl-mode', 'raw', '--scope', 'combined_chromosomes',
        '--output-prefix', str(prefix),
    ]) == 0
    category_roots = list(tmp_path.glob('*_scorepct0_states'))
    assert len(category_roots) == 1
    root = category_roots[0]
    expected = {
        'active_genes': ({1: (186, 6), 2: (372, 4), 3: (558, 2)}, 186),
        'repressed_genes': ({1: (191, 4), 2: (382, 3), 3: (573, 2)}, 191),
    }
    for state, (counts_by_order, expected_nrl) in expected.items():
        folder = root / state
        assert (folder / 'distance_distribution.png').is_file()
        assert (folder / 'summary.tsv').is_file()
        summary = list(csv.DictReader((folder / 'summary.tsv').open(), delimiter='\t'))
        assert {int(row['order']) for row in summary} == {1, 2, 3}
        for row in summary:
            order = int(row['order'])
            expected_distance, expected_count = counts_by_order[order]
            assert row['state'] == state
            assert int(row['total_pairs']) == expected_count
            assert int(row['raw_mode_bp']) == expected_distance
        regressions = list(csv.DictReader(
            (folder / 'state_nrl_regression_summary.tsv').open(), delimiter='\t'))
        assert len(regressions) == 1
        assert regressions[0]['state'] == state
        assert float(regressions[0]['nrl_bp']) == pytest.approx(expected_nrl)
        assert (folder / 'state_nrl_regression_combined_chromosomes.png').is_file()
    # The pooled +1 set still includes pairs crossing interval boundaries.
    pooled = next(tmp_path.glob('*_scorepct0_summary.tsv'))
    rows = list(csv.DictReader(pooled.open(), delimiter='\t'))
    all_one = next(row for row in rows if row['state'] == 'All' and row['order'] == '1')
    assert int(all_one['total_pairs']) == 12
    # Different intervals of the SAME label may not form a +1 pair.
    active_one = next(row for row in rows if row['state'] == 'active_genes' and row['order'] == '1')
    assert int(active_one['total_pairs']) == 6


def test_overlapping_same_category_intervals_do_not_double_count_pairs(tmp_path):
    peaks, states = _write_example(tmp_path, overlapping=True)
    indexes, _ = distances.build_state_indexes(states)
    by_chrom, _, _ = distances.load_peaks(peaks, state_indexes=indexes, position_column=7)
    chrom, genome = distances.compute_within_interval_state_counts(
        by_chrom, indexes, max_order=3,
    )
    assert genome[1]['active_genes'] == Counter({186: 6})
    assert genome[2]['active_genes'] == Counter({372: 4})
    assert genome[3]['active_genes'] == Counter({558: 2})
    assert chrom[1]['1']['active_genes'] == Counter({186: 6})


def test_category_order_restarts_and_contig_aliases(tmp_path):
    peaks, states = _write_example(tmp_path)
    indexes, _ = distances.build_state_indexes(states)
    by_chrom, _, _ = distances.load_peaks(peaks, state_indexes=indexes, position_column=7)
    chrom, genome = distances.compute_within_interval_state_counts(
        by_chrom, indexes, max_order=8,
    )
    assert '1' in chrom[1]
    assert 4 not in genome or 'active_genes' not in genome[4]
    # Repressed has five sites, so the fourth-order distance remains within it.
    assert genome[4]['repressed_genes'] == Counter({764: 1})


def test_legacy_endpoint_selection_remains_explicit(tmp_path):
    peaks, states = _write_example(tmp_path)
    indexes, _ = distances.build_state_indexes(states)
    by_chrom, _, _ = distances.load_peaks(peaks, state_indexes=indexes, position_column=7)
    pooled = distances.compute_distance_counts(
        by_chrom, threshold=0, min_distance=1, max_distance=1000,
        max_order=1, duplicate_policy='highest-score',
    )
    assert pooled.genome_state[1]['active_genes'][242] == 1
