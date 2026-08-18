from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import numpy as np

from nucleosuite.core.regions import ProcessingRegion
from nucleosuite.scoring import basic_tracks
from nucleosuite.workflows import tracks


class FakeContext:
    def __init__(self, fragments):
        self.fragments = list(fragments)
        self.regions = [ProcessingRegion("chr1", 0, 300, 0, 300)]
        self.bigwig_header = [("chr1", 300)]
        self.collect_calls = 0
        self.closed = False

    def collect(self, **_kwargs):
        self.collect_calls += 1
        return list(self.fragments)

    def close(self):
        self.closed = True


def make_args(tmp_path: Path, specs: list[str], *, cap: int = 0):
    return Namespace(
        bamfiles=["sample.bam"],
        fragment_files=None,
        fasta=None,
        chrom_sizes=None,
        contigs=["chr1"],
        fragment_range=specs,
        spec_file=[],
        output_dir=str(tmp_path),
        output_prefix="sample",
        report=str(tmp_path / "report.tsv"),
        max_duplicates=1,
        max_per_coordinate=cap,
        dedup_scope="all_bams",
        chunk_bp=100_000,
        overlap_bp=1_000,
        subsample=None,
        seed=1,
        even_dyad="split",
        output_format="none",
        interval_format="bed",
        pns_mode_length=167,
        pns_smooth_window=0,
        pns_smooth_order=2,
        pns_min_region_length=50,
        pns_max_neg_run=0,
        pns_peak_score_scale=1.0,
        wps_protection=120,
        wps_baseline_window=1000,
        wps_sg_window=21,
        wps_sg_order=2,
        wps_peak_track="sm_mWPS",
        wps_peak_score_scale=1.0,
        wps_peak_minlen=50,
        wps_peak_maxlen=150,
        wps_peak_maxregion=450,
        wps_peak_merge_gap=5,
        wps_peak_varicutoff=5.0,
    )


def test_parse_single_and_bounded_fragment_ranges():
    assert tracks.parse_fragment_range("145") == tracks.FragmentRange(145, 145)
    assert tracks.parse_fragment_range("137-197") == tracks.FragmentRange(137, 197)


def test_combined_engine_fetches_once_and_assigns_overlapping_ranges(monkeypatch, tmp_path):
    context = FakeContext([(10, 155), (20, 181)])  # 145 bp and 161 bp
    monkeypatch.setattr(tracks, "prepare_fragment_run", lambda **_kwargs: context)
    monkeypatch.setattr(tracks, "open_track_handles", lambda **_kwargs: ({}, {}))
    monkeypatch.setattr(tracks, "close_track_handles", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(tracks, "write_fragment_outputs", lambda **_kwargs: ("", ""))

    captured = {}

    def capture(**kwargs):
        prefix = next(
            spec.output_prefix
            for spec in active_specs
            if spec.output_tracks == kwargs["tracks"]
            and spec.output_prefix not in captured
        )
        captured[prefix] = {
            name: np.asarray(kwargs["scores"][name][0][2]).copy()
            for name in kwargs["tracks"]
        }

    monkeypatch.setattr(tracks, "write_tracks", capture)
    args = make_args(
        tmp_path,
        [
            "137-197=coverage,dyad,fragment_left_ends,fragment_right_ends",
            "145=dyad",
        ],
    )
    active_specs = tracks.load_specs(args)
    monkeypatch.setattr(tracks, "load_specs", lambda _args: active_specs)

    assert tracks.run(args) == 0
    assert context.collect_calls == 1
    broad = captured[active_specs[0].output_prefix]
    exact = captured[active_specs[1].output_prefix]
    assert broad["coverage"].sum() == 145 + 161
    assert broad["fragment_left_ends"].sum() == 2
    assert broad["fragment_right_ends"].sum() == 2
    assert broad["dyad"].sum() == 2
    assert exact["dyad"].sum() == 1
    assert Path(args.report).is_file()


def test_sparse_coordinate_cap_defaults_to_unlimited():
    arrays = basic_tracks.new_arrays(300)
    basic_tracks.add_fragment(arrays, 10, 155, 0, 300)
    basic_tracks.add_fragment(arrays, 10, 156, 0, 300)
    assert arrays["fragment_left_ends"][10] == 2
    basic_tracks.cap_sparse_arrays(arrays, 0)
    assert arrays["fragment_left_ends"][10] == 2
    basic_tracks.cap_sparse_arrays(arrays, 1)
    assert arrays["fragment_left_ends"][10] == 1


def test_spec_file_can_use_explicit_output_prefixes(tmp_path):
    spec_file = tmp_path / "spec.tsv"
    prefix = tmp_path / "01_pns" / "sample_PNS_mode167_lower137_upper197"
    spec_file.write_text(
        "fragment_range\toutput_prefix\ttracks\tbasic_scope\n"
        f"137-197\t{prefix}\tpns,posPNS,coverage,pns_peaks\trange\n"
    )
    args = make_args(tmp_path, [])
    args.spec_file = [str(spec_file)]
    specs = tracks.load_specs(args)
    assert len(specs) == 1
    assert specs[0].output_prefix == str(prefix)
    assert specs[0].basic_scope == "range"
    assert specs[0].tracks == ("pns", "posPNS", "coverage", "pns_peaks")


def test_combined_pns_and_wps_match_direct_scoring(monkeypatch, tmp_path):
    from nucleosuite.scoring import pns as pns_scoring
    from nucleosuite.scoring import wps as wps_scoring

    fragments = [(40, 185), (80, 241)]  # 145 bp and 161 bp
    context = FakeContext(fragments)
    monkeypatch.setattr(tracks, "prepare_fragment_run", lambda **_kwargs: context)
    monkeypatch.setattr(tracks, "open_track_handles", lambda **_kwargs: ({}, {}))
    monkeypatch.setattr(tracks, "close_track_handles", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(tracks, "write_fragment_outputs", lambda **_kwargs: ("", ""))
    captured = {}

    def capture(**kwargs):
        for name in kwargs["tracks"]:
            captured[name] = np.asarray(kwargs["scores"][name][0][2]).copy()

    monkeypatch.setattr(tracks, "write_tracks", capture)
    args = make_args(tmp_path, ["137-197=pns,posPNS,wps,coverage"])
    args.wps_baseline_window = 100
    active_specs = tracks.load_specs(args)
    monkeypatch.setattr(tracks, "load_specs", lambda _args: active_specs)
    tracks.run(args)

    basic = basic_tracks.new_arrays(300)
    pns_arrays = pns_scoring.new_arrays(300)
    centred, positive = pns_scoring.precompute_distributions(range(137, 198), 167)
    wps_array = wps_scoring.new_array(300)
    wps_distributions = wps_scoring.precompute_distributions(range(137, 198), 120)
    for start, end in fragments:
        basic_tracks.add_fragment(basic, start, end, 0, 300, "split")
        pns_scoring.add_fragment(
            pns_arrays, start, end, 0, 300, 167, centred, positive
        )
        wps_scoring.add_fragment(
            wps_array, start, end, 0, 120, wps_distributions
        )
    expected_pns = pns_scoring.to_scores(pns_arrays, "chr1", 0, 0, 2)
    expected_wps = wps_scoring.to_scores(wps_array, "chr1", 0, 100, 21, 2)
    assert np.array_equal(captured["coverage"], basic["coverage"])
    assert np.allclose(captured["pns"], expected_pns["pns"][0][2])
    assert np.allclose(captured["posPNS"], expected_pns["posPNS"][0][2])
    assert np.allclose(captured["wps"], expected_wps["wps"][0][2])


def test_tracks_cli_defaults_leave_sparse_coordinate_values_unlimited():
    from nucleosuite.cli.main import build_parser

    args = build_parser().parse_args([
        "tracks", "--bam", "sample.bam",
        "--fragment-range", "145=dyad",
    ])
    assert args.max_duplicates == 1
    assert args.max_per_coordinate == 0


def test_wps_peak_overlap_is_validated_after_specs_load(tmp_path):
    args = make_args(tmp_path, ["120-180=wps_peaks"])
    args.overlap_bp = 959
    try:
        tracks.run(args)
    except ValueError as error:
        assert "--overlap-bp >= 960" in str(error)
    else:
        raise AssertionError("Expected WPS overlap validation to fail")


def test_combined_sequence_features_are_computed_once_per_fragment(monkeypatch, tmp_path):
    context = FakeContext([(10, 155)])
    context.fasta = object()
    spec_file = tmp_path / "sequence_specs.tsv"
    spec_file.write_text(
        "fragment_range\toutput_prefix\ttracks\tbasic_scope\n"
        f"145\t{tmp_path / 'dinuc'}\tdinuc_profile\trange\n"
        f"145\t{tmp_path / 'types'}\tww_types\trange\n"
        f"145\t{tmp_path / 'type_dyads'}\ttype_dyads\trange\n"
    )
    args = make_args(tmp_path, [])
    args.fasta = "genome.fa"
    args.spec_file = [str(spec_file)]

    monkeypatch.setattr(tracks, "prepare_fragment_run", lambda **_kwargs: context)
    monkeypatch.setattr(tracks, "prepare_reference_if_needed", lambda **_kwargs: object())
    monkeypatch.setattr(tracks, "open_track_handles", lambda **_kwargs: ({}, {}))
    monkeypatch.setattr(tracks, "close_track_handles", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(tracks, "write_tracks", lambda **_kwargs: None)
    monkeypatch.setattr(tracks, "write_fragment_outputs", lambda **_kwargs: ("", ""))
    monkeypatch.setattr(tracks, "write_fragment_bed_rows", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(tracks, "_write_sequence_outputs", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(tracks, "finalise_interval_files", lambda *_args, **_kwargs: [])

    calls = []
    def features(**kwargs):
        calls.append((kwargs["fragment_start"], kwargs["fragment_end"]))
        return "A" * 145, "type1"
    monkeypatch.setattr(tracks, "_sequence_features_for_fragment", features)

    assert tracks.run(args) == 0
    assert calls == [(10, 155)]


def test_combined_sequence_feature_path_uses_real_core_classifier(monkeypatch):
    monkeypatch.setattr(
        tracks,
        "extract_reference_sequence",
        lambda **kwargs: "A" * (kwargs["seq_end"] - kwargs["seq_start"]),
    )
    fragment_sequence, ww_type = tracks._sequence_features_for_fragment(
        fasta=object(),
        reference_context=object(),
        fragment_start=10,
        fragment_end=155,
        need_ww_type=True,
    )
    assert fragment_sequence == "A" * 145
    assert ww_type == "type1"


def test_combined_ww_type_beds_are_coordinate_sorted(monkeypatch, tmp_path):
    context = FakeContext([(20, 165), (10, 155), (10, 154)])
    context.fasta = object()
    prefix = tmp_path / "ww" / "sample_wwtypes_lower144_upper145"
    spec_file = tmp_path / "sequence_specs.tsv"
    spec_file.write_text(
        "fragment_range\toutput_prefix\ttracks\tbasic_scope\n"
        f"144-145\t{prefix}\tww_types\trange\n",
        encoding="utf-8",
    )
    args = make_args(tmp_path, [])
    args.fasta = "genome.fa"
    args.spec_file = [str(spec_file)]

    monkeypatch.setattr(tracks, "prepare_fragment_run", lambda **_kwargs: context)
    monkeypatch.setattr(tracks, "prepare_reference_if_needed", lambda **_kwargs: object())
    monkeypatch.setattr(tracks, "open_track_handles", lambda **_kwargs: ({}, {}))
    monkeypatch.setattr(tracks, "close_track_handles", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(tracks, "write_tracks", lambda **_kwargs: None)
    monkeypatch.setattr(tracks, "write_fragment_outputs", lambda **_kwargs: ("", ""))
    monkeypatch.setattr(tracks, "_write_sequence_outputs", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(tracks, "finalise_interval_files", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(
        tracks,
        "_sequence_features_for_fragment",
        lambda **kwargs: (
            "A" * (kwargs["fragment_end"] - kwargs["fragment_start"]),
            "type1",
        ),
    )

    assert tracks.run(args) == 0
    assert Path(f"{prefix}_ww_types.bed").read_text(encoding="utf-8").splitlines() == [
        "chr1\t10\t154\ttype1",
        "chr1\t10\t155\ttype1",
        "chr1\t20\t165\ttype1",
    ]
    assert Path(f"{prefix}_type1.bed").read_text(encoding="utf-8").splitlines() == [
        "chr1\t10\t154",
        "chr1\t10\t155",
        "chr1\t20\t165",
    ]


def test_staged_bedgraph_writer_validates_during_generation(tmp_path):
    from nucleosuite.io.bedgraph import ValidatedBedGraphWriter
    import json

    path = tmp_path / "staged" / "coverage.bedGraph"
    writer = ValidatedBedGraphWriter(
        path,
        track="coverage",
        chrom_order=["chr1"],
        source_bigwig=tmp_path / "coverage.bw",
    )
    writer.add_values(
        chrom="chr1",
        start=0,
        values=np.asarray([0.0, 0.0, 1.0, 1.0, 0.0]),
        sparse=False,
    )
    writer.add_values(
        chrom="chr1",
        start=5,
        values=np.asarray([0.0, 2.0]),
        sparse=False,
    )
    writer.close()

    assert path.read_text(encoding="utf-8").splitlines() == [
        "chr1\t0\t2\t0",
        "chr1\t2\t4\t1",
        "chr1\t4\t6\t0",
        "chr1\t6\t7\t2",
    ]
    metadata = json.loads(
        path.with_name(path.name + ".complete.json").read_text(encoding="utf-8")
    )
    assert metadata["complete"] is True
    assert metadata["sorted"] is True
    assert metadata["nonoverlapping"] is True
    assert metadata["records"] == 4


def test_staged_bedgraph_writer_rejects_coordinate_overlap(tmp_path):
    from nucleosuite.io.bedgraph import ValidatedBedGraphWriter
    import pytest

    path = tmp_path / "overlap.bedGraph"
    writer = ValidatedBedGraphWriter(
        path,
        track="coverage",
        chrom_order=["chr1"],
        source_bigwig=tmp_path / "coverage.bw",
    )
    writer.add_interval("chr1", 10, 20, 1.0)
    writer.add_interval("chr1", 15, 25, 2.0)
    with pytest.raises(ValueError, match="Overlapping bedGraph intervals"):
        writer.close()
    assert not path.exists()
    assert not path.with_name(path.name + ".complete.json").exists()


def test_tracks_workflow_stages_bedgraph_under_combined_root(monkeypatch, tmp_path):
    context = FakeContext([(10, 15)])
    worker_root = tmp_path / "per_contig" / "chr1"
    combined_root = tmp_path / "combined"
    args = make_args(worker_root, ["5=coverage"])
    args.output_dir = str(worker_root / "01_combined_tracks")
    args.staged_bedgraph_root = str(
        combined_root / "temporary_bedgraph_combine"
    )
    args.staged_bedgraph_source_root = str(worker_root)
    args.staged_bedgraph_source_id = "chr1"

    monkeypatch.setattr(tracks, "prepare_fragment_run", lambda **_kwargs: context)
    monkeypatch.setattr(tracks, "open_track_handles", lambda **_kwargs: ({}, {}))
    monkeypatch.setattr(tracks, "close_track_handles", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(tracks, "write_fragment_outputs", lambda **_kwargs: ("", ""))

    assert tracks.run(args) == 0
    staged = (
        combined_root
        / "temporary_bedgraph_combine"
        / "per_contig"
        / "chr1"
        / "01_combined_tracks"
        / "exact_5"
        / "sample_exact_5_coverage.bedGraph"
    )
    assert staged.is_file()
    assert staged.read_text(encoding="utf-8").splitlines() == [
        "chr1\t0\t10\t0",
        "chr1\t10\t15\t1",
        "chr1\t15\t300\t0",
    ]
    assert staged.with_name(staged.name + ".complete.json").is_file()
