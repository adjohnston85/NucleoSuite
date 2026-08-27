import numpy as np

from nucleosuite.peaks.pns import (
    call_records,
    filter_records_by_coverage,
    find_peaks_and_regions,
    write_records,
)
from nucleosuite.peaks.wps import call_records as call_wps_records


def test_pns_peak_core_ownership():
    signal = np.r_[np.zeros(5), np.ones(60), np.zeros(5)]
    records = call_records(signal, "ctg", 100, 100, 170, 50, 5)
    assert len(records) == 1
    assert records[0]["chrom"] == "ctg"


def test_wps_preserves_contig_name():
    signal = np.r_[np.zeros(5), np.full(60, 10.0), np.zeros(5)]
    records = call_wps_records(
        signal, "scaffold_1", 0, 0, len(signal),
        min_length=50, max_length=150, score_cutoff=5,
    )
    assert records
    assert records[0]["chrom"] == "scaffold_1"


def test_peakcall_cli_aliases_and_option_names():
    from nucleosuite.cli.main import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "peakcall", "--input-bigwig", "signal.bw",
        "--peak-caller", "wps", "--call-type", "nucleosome",
    ])
    assert args.method == "wps"
    assert args.signal == "nucleosome"


def test_wps_peak_caller_is_default():
    from nucleosuite.cli.main import build_parser

    parser = build_parser()
    args = parser.parse_args(["wps", "-b", "sample.bam"])
    assert args.peak_caller == "wps"
    assert args.protection == "auto"
    assert args.mode_histogram_smoothing == "none"


def test_call_peaks_method_specific_smoothing_defaults():
    from nucleosuite.cli.call_peaks import resolve_smooth_window

    assert resolve_smooth_window("pns", "raw", None) == 0
    assert resolve_smooth_window("wps", "raw", None) == 21
    assert resolve_smooth_window("wps", "adjusted", None) == 0
    assert resolve_smooth_window("wps", "raw", 0) == 0


def test_pns_peak_defaults_do_not_bridge_nonpositive_signal():
    signal = np.r_[np.ones(30), np.zeros(1), np.ones(30)]
    call = find_peaks_and_regions(signal, genomic_start=0, min_length=20)
    assert call["regions"] == [(0, 30), (31, 61)]


def test_nuc_score_cli_defaults_to_sns_raw_tracks_and_zero_gap():
    from nucleosuite.cli.main import build_parser

    args = build_parser().parse_args(["nuc-score", "--bam", "sample.bam"])
    assert args.mode_length == "auto"
    assert args.mode_histogram_smoothing == "none"
    assert args.smooth_window == 0
    assert args.max_neg_run == 0
    assert args.scoring_method == "sns"
    assert args.score_tracks == ["sns", "posSNS"]
    assert args.bigbed_score_scale == 1000.0




def test_nuc_score_cli_exposes_pns_as_optional_scoring_method():
    from nucleosuite.cli.main import build_parser

    args = build_parser().parse_args([
        "nuc-score", "--bam", "sample.bam", "--scoring-method", "pns", "--mode", "167"
    ])
    assert args.scoring_method == "pns"
    assert args.score_tracks == ["sns", "posSNS"]


def test_nuc_score_default_sns_run_keeps_sns_tracks(monkeypatch, tmp_path):
    from nucleosuite.cli.main import build_parser
    from nucleosuite.cli import nuc_score as nuc_score_cli
    import nucleosuite.parallel as parallel

    args = build_parser().parse_args([
        "nuc-score", "--bam", "sample.bam", "--mode", "167",
        "--out-prefix", str(tmp_path / "sample"),
    ])
    captured = {}

    def fake_run_native(command, run_args, runner):
        captured["command"] = command
        captured["method"] = run_args.scoring_method
        captured["tracks"] = list(run_args.score_tracks)
        return 0

    monkeypatch.setattr(parallel, "run_native_per_contig", fake_run_native)
    assert nuc_score_cli.run(args) == 0
    assert captured == {
        "command": "nuc-score",
        "method": "sns",
        "tracks": ["sns", "posSNS"],
    }


def test_nuc_score_cli_exposes_bns_as_optional_scoring_method():
    from nucleosuite.cli.main import build_parser

    args = build_parser().parse_args([
        "nuc-score", "--bam", "sample.bam", "--scoring-method", "bns"
    ])
    assert args.scoring_method == "bns"
    assert args.score_tracks == ["sns", "posSNS"]


def test_nuc_score_cli_exposes_tns_as_optional_scoring_method():
    from nucleosuite.cli.main import build_parser

    args = build_parser().parse_args([
        "nuc-score", "--bam", "sample.bam", "--scoring-method", "tns", "--mode", "167"
    ])
    assert args.scoring_method == "tns"
    assert args.score_tracks == ["sns", "posSNS"]


def test_nuc_score_tns_run_maps_default_score_tracks(monkeypatch, tmp_path):
    from nucleosuite.cli.main import build_parser
    from nucleosuite.cli import nuc_score as nuc_score_cli
    import nucleosuite.parallel as parallel

    args = build_parser().parse_args([
        "nuc-score", "--bam", "sample.bam", "--scoring-method", "tns", "--mode", "167",
        "--out-prefix", str(tmp_path / "sample"),
    ])
    captured = {}

    def fake_run_native(command, run_args, runner):
        captured["command"] = command
        captured["tracks"] = list(run_args.score_tracks)
        return 0

    monkeypatch.setattr(parallel, "run_native_per_contig", fake_run_native)
    assert nuc_score_cli.run(args) == 0
    assert captured["command"] == "nuc-score"
    assert captured["tracks"] == ["tns", "posTNS"]


def test_nuc_score_cli_peak_coverage_filter_is_off_by_default_and_configurable():
    from nucleosuite.cli.main import build_parser

    parser = build_parser()
    default = parser.parse_args(["nuc-score", "--bam", "sample.bam"])
    selected = parser.parse_args([
        "nuc-score", "--bam", "sample.bam", "--peak-coverage-threshold", "2"
    ])
    assert default.peak_coverage_threshold is None
    assert selected.peak_coverage_threshold == 2.0


def test_nuc_score_peak_calling_can_be_disabled_without_disabling_score_tracks():
    from nucleosuite.cli.main import build_parser

    parser = build_parser()
    default = parser.parse_args(["nuc-score", "--bam", "sample.bam"])
    tracks_only = parser.parse_args(
        ["nuc-score", "--bam", "sample.bam", "--no-peak-calling"]
    )
    assert default.peak_calling is True
    assert tracks_only.peak_calling is False
    assert tracks_only.score_mode == "on"


def test_pns_peak_coverage_filter_uses_bed_column7_position():
    coverage = np.zeros(30, dtype=float)
    coverage[15] = 2.0
    coverage[20] = 1.0
    records = [
        {"region_centre": 115, "chrom": "chr1"},
        {"region_centre": 120, "chrom": "chr1"},
    ]
    retained, filtered = filter_records_by_coverage(
        records, coverage, coverage_start=100, threshold=2.0
    )
    assert retained == [records[0]]
    assert filtered == 1


def test_nuc_score_bigbed_score_scale_is_configurable():
    from nucleosuite.cli.main import build_parser

    parser = build_parser()
    pns = parser.parse_args([
        "nuc-score", "--bam", "sample.bam", "--bigbed-score-scale", "250"
    ])
    called = parser.parse_args([
        "call-peaks", "--input-bigwig", "signal.bw", "--peak-caller", "pns",
        "--bigbed-score-scale", "500",
    ])
    tracks = parser.parse_args([
        "tracks", "--bam", "sample.bam", "--fragment-range", "137-197=pns_peaks",
        "--bigbed-score-scale", "750",
    ])
    assert pns.bigbed_score_scale == 250.0
    assert called.bigbed_score_scale == 500.0
    assert tracks.bigbed_score_scale == 750.0


def test_pns_bed_scores_keep_six_decimal_places(tmp_path):
    path = tmp_path / "pns_peaks.bed"
    write_records(
        str(path),
        [
            {
                "chrom": "chr1",
                "region_start": 10,
                "region_end": 20,
                "region_centre": 15,
                "peak_score": 3.14159265,
            }
        ],
        "nuc",
        1.0,
        "w",
    )
    assert path.read_text().split("\t")[4] == "3.141593"
