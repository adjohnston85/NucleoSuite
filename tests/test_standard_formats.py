"""Tests for standard public output formats."""

from __future__ import annotations

from collections import Counter

from nucleosuite.distances import (
    DistanceResults,
    ThresholdSelection,
    PeakRecord,
    build_parser as build_distances_parser,
    compute_distance_counts,
    load_peaks,
    write_distribution_outputs,
    write_filtered_bed,
    write_threshold_metadata,
)
from nucleosuite.io.summaries import print_fragment_length_histogram, write_fragment_outputs
from nucleosuite.peaks.pns import write_records


def test_fragment_summary_is_headered_tsv(tmp_path):
    prefix = str(tmp_path / "sample")
    summary_path, lengths_path = write_fragment_outputs(
        prefix,
        total_fragments_filtered=10,
        total_fragments_used=8,
        unique_bases_covered=500,
        length_counts=Counter({147: 5, 167: 3}),
        dedup_scope="all_bams",
        max_duplicates=1,
        max_per_coordinate=0,
    )
    assert summary_path.endswith("_fragment_summary.tsv")
    assert lengths_path.endswith("_fragment_length_counts.tsv")
    assert (tmp_path / "sample_fragment_summary.tsv").read_text().splitlines()[0] == "metric\tvalue"
    assert (tmp_path / "sample_fragment_length_counts.tsv").read_text().splitlines()[0] == "fragment_length\tcount"
    assert (tmp_path / "sample_fragment_length_distribution.png").exists()


def test_pns_peak_writer_outputs_bed8(tmp_path):
    output = tmp_path / "peaks.bed"
    records = [{
        "chrom": "chr1",
        "region_start": 100,
        "region_end": 180,
        "region_centre": 140,
        "peak_score": 12.4,
    }]
    write_records(str(output), records, "nuc", score_scale=10.0, mode="w")
    fields = output.read_text().rstrip().split("\t")
    assert len(fields) == 8
    assert fields[:3] == ["chr1", "100", "180"]
    assert fields[4] == "124.000000"
    assert fields[5:] == [".", "140", "141"]


def test_distances_uses_bed_score_column_and_writes_bed8(tmp_path):
    input_bed = tmp_path / "input.bed"
    input_bed.write_text(
        "chr1\t100\t180\tpeak1\t250\t.\t140\t141\n"
        "chr1\t300\t380\tpeak2\t500\t+\t340\t341\n"
    )
    peaks, scores, summary = load_peaks(input_bed, state_indexes=None)
    assert scores.tolist() == [250.0, 500.0]
    assert summary.used_lines == 2

    output_bed = tmp_path / "filtered.bed"
    assert write_filtered_bed(peaks, output_bed) == 2
    lines = output_bed.read_text().splitlines()
    assert all(len(line.split("\t")) == 8 for line in lines)
    assert lines[0].split("\t")[6:8] == ["140", "141"]


def test_distances_state_bed_is_option_only():
    parser = build_distances_parser()
    args = parser.parse_args(["peaks.bed", "--state-bed", "states.bed"])
    assert args.score_column == 5
    assert args.state_bed == "states.bed"


def test_distance_tables_are_headered_tsv(tmp_path):
    results = DistanceResults(
        chrom_state={},
        chrom_all={1: {"chr1": Counter({180: 2, 181: 1})}},
        genome_state={},
        genome_all={1: Counter({180: 2, 181: 1})},
        duplicates={},
        retained_by_chrom={},
        threshold_pass_count=4,
        retained_count=4,
    )
    distance_path = tmp_path / "distances.tsv"
    summary_path = tmp_path / "summary.tsv"
    metadata_path = tmp_path / "metadata.tsv"
    selection = ThresholdSelection(50.0, 50.0, 100.0)

    write_threshold_metadata(
        metadata_path,
        input_path="peaks.bed",
        state_path=None,
        selection=selection,
        results=results,
        score_column=5,
        min_distance=1,
        max_distance=1000,
        max_order=1,
        duplicate_policy="highest-score",
    )
    write_distribution_outputs(
        results,
        distance_path=distance_path,
        summary_path=summary_path,
        max_order=1,
        include_chromosomes=True,
        include_genome=True,
        include_state_strata=False,
        include_zero_distances=False,
        count_smooth_window=0,
        count_smooth_polyorder=2,
        percent_smooth_window=0,
        percent_smooth_polyorder=3,
    )

    assert metadata_path.read_text().splitlines()[0] == "parameter\tvalue"
    assert distance_path.read_text().splitlines()[0].startswith("order\tscope")
    assert summary_path.read_text().splitlines()[0].startswith("order\tscope")


def test_distance_counting_accepts_bed8_peak_records():
    peaks = {
        "chr1": [
            PeakRecord(140, 250.0, "Active", 100, 180, "peak1", "."),
            PeakRecord(340, 500.0, "Active", 300, 380, "peak2", "+"),
            PeakRecord(550, 300.0, "Repressed", 500, 600, "peak3", "-"),
        ]
    }
    results = compute_distance_counts(
        peaks, threshold=0.0, min_distance=1, max_distance=1000,
        max_order=2, duplicate_policy="highest-score",
    )
    assert results.genome_all[1][200] == 1
    assert results.genome_all[1][210] == 1
    assert results.genome_all[2][410] == 1
    assert results.genome_state[1]["Active"][200] == 1


def test_console_fragment_histogram(capsys):
    print_fragment_length_histogram(Counter({145: 5, 147: 10}), label="sample", width=10)
    output = capsys.readouterr().out
    assert "Fragment-length distribution: sample" in output
    assert "145" in output and "-----" in output
    assert "147" in output and "----------" in output
