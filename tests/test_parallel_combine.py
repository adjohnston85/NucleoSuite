from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pytest

from nucleosuite import combine
from nucleosuite.dac import write_dac_tsv
from nucleosuite.sequence.dinucleotide import DINUCS


def test_combine_resource_defaults_keep_indexed_writers_independent() -> None:
    args = combine.build_parser().parse_args(["--input-dir", "run", "--cores", "8"])
    assert args.cores == 8
    assert args.streaming_combine_cores is None
    assert args.indexed_combine_cores == 1
    assert args.combine_chunk_bp == 100_000

    overridden = combine.build_parser().parse_args(
        [
            "--input-dir", "run",
            "--cores", "8",
            "--streaming-combine-cores", "3",
            "--indexed-combine-cores", "2",
            "--combine-chunk-bp", "4096",
        ]
    )
    assert overridden.streaming_combine_cores == 3
    assert overridden.indexed_combine_cores == 2
    assert overridden.combine_chunk_bp == 4096


def _write_dac_summary(path: Path, output: Path, total_signal: float) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["State", "Output", "Total signal"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"State": "all", "Output": str(output), "Total signal": total_signal})


def test_dac_combination_sums_raw_values_and_opportunities(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(combine, "plot_dac_tsv", lambda *args, **kwargs: None)
    input_a = tmp_path / "chr1" / "sample_chr1_DAC_opportunity_normalized.tsv"
    input_b = tmp_path / "chr2" / "sample_chr2_DAC_opportunity_normalized.tsv"
    input_a.parent.mkdir()
    input_b.parent.mkdir()
    write_dac_tsv(
        str(input_a),
        np.array([0.0, 10.0, 20.0]),
        np.array([0.0, 2.0, 4.0]),
        True,
        10.0,
        1_000_000.0,
    )
    write_dac_tsv(
        str(input_b),
        np.array([0.0, 5.0, 30.0]),
        np.array([0.0, 3.0, 6.0]),
        True,
        20.0,
        1_000_000.0,
    )
    _write_dac_summary(input_a.parent / "sample_chr1_DAC_summary.tsv", input_a, 10.0)
    _write_dac_summary(input_b.parent / "sample_chr2_DAC_summary.tsv", input_b, 20.0)

    output = tmp_path / "combined" / "sample_DAC_opportunity_normalized.tsv"
    combine._combine_dac([input_a, input_b], output)

    fields, rows = combine._read_tsv(output)
    assert fields[:5] == [
        "Distance",
        "DAC Value",
        "DAC Value Percent",
        "Raw DAC Value",
        "Opportunities",
    ]
    by_distance = {int(row["Distance"]): row for row in rows}
    assert float(by_distance[1]["Raw DAC Value"]) == pytest.approx(15.0)
    assert float(by_distance[1]["Opportunities"]) == pytest.approx(5.0)
    assert float(by_distance[1]["DAC Value"]) == pytest.approx(3.0)
    assert float(by_distance[2]["Raw DAC Value"]) == pytest.approx(50.0)
    assert float(by_distance[2]["Opportunities"]) == pytest.approx(10.0)
    assert float(by_distance[2]["DAC Value"]) == pytest.approx(5.0)
    assert float(by_distance[1]["DAC Value Percent"]) == pytest.approx(37.5)
    assert float(by_distance[2]["DAC Value Percent"]) == pytest.approx(62.5)
    assert float(by_distance[1]["DAC per million signal-pairs"]) == pytest.approx(
        15.0 / (30.0 * 30.0) * 1_000_000.0
    )


def _write_dinuc_counts(path: Path, n_valid: int, aa: int, tt: int, used: int) -> None:
    fields = ["position", "n_valid", *[f"{name}_count" for name in DINUCS], "fragments_used", "fragments_skipped"]
    row = {name: 0 for name in fields}
    row.update(
        {
            "position": 0,
            "n_valid": n_valid,
            "AA_count": aa,
            "TT_count": tt,
            "fragments_used": used,
            "fragments_skipped": 1,
        }
    )
    combine._write_tsv(path, fields, [row])


def test_dinucleotide_combination_uses_counts_and_valid_positions(tmp_path: Path) -> None:
    first = tmp_path / "a_counts.tsv"
    second = tmp_path / "b_counts.tsv"
    output = tmp_path / "combined_counts.tsv"
    profile = tmp_path / "combined.tsv"
    _write_dinuc_counts(first, n_valid=10, aa=4, tt=1, used=3)
    _write_dinuc_counts(second, n_valid=30, aa=6, tt=9, used=7)

    combine._combine_dinuc_counts([first, second], output)
    combine._profile_from_counts(output, profile, fraction=False)

    _, rows = combine._read_tsv(output)
    assert len(rows) == 1
    assert int(rows[0]["n_valid"]) == 40
    assert int(rows[0]["AA_count"]) == 10
    assert int(rows[0]["TT_count"]) == 10
    assert int(rows[0]["fragments_used"]) == 10
    assert int(rows[0]["fragments_skipped"]) == 2

    _, profile_rows = combine._read_tsv(profile)
    assert float(profile_rows[0]["AA_pct"]) == pytest.approx(25.0)
    assert float(profile_rows[0]["TT_pct"]) == pytest.approx(25.0)
    assert float(profile_rows[0]["WW_pct"]) == pytest.approx(50.0)


def _write_converter(path: Path, *, fail: bool) -> None:
    body = "#!/usr/bin/env bash\n"
    if fail:
        body += "exit 7\n"
    else:
        body += 'cp "$1" "$3"\n'
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)


def test_bigwig_temp_bedgraphs_removed_only_after_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    converter = bin_dir / "bedGraphToBigWig"
    _write_converter(converter, fail=False)
    monkeypatch.setenv("PATH", f"{bin_dir}:{Path('/usr/bin')}")

    def fake_bedgraph(_input: Path, output: Path, _order, **_kwargs) -> tuple[int, int]:
        output.write_text("chr1\t0\t2\t1\n", encoding="utf-8")
        return 1, 1

    monkeypatch.setattr(combine, "_bigwig_to_bedgraph", fake_bedgraph)
    monkeypatch.setattr(combine, "_verify_bigwig", lambda *args, **kwargs: None)
    temp_dir = tmp_path / "tmp_bedgraphs"
    output = tmp_path / "combined.bw"
    combine._combine_bigwig(
        [tmp_path / "chr1.bw", tmp_path / "chr2.bw"],
        output,
        [("chr1", 10), ("chr2", 10)],
        temp_dir,
    )
    assert output.exists()
    assert not temp_dir.exists()


def test_bigwig_temp_bedgraphs_retained_after_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    converter = bin_dir / "bedGraphToBigWig"
    _write_converter(converter, fail=True)
    monkeypatch.setenv("PATH", f"{bin_dir}:{Path('/usr/bin')}")

    def fake_bedgraph(_input: Path, output: Path, _order, **_kwargs) -> tuple[int, int]:
        output.write_text("chr1\t0\t2\t1\n", encoding="utf-8")
        return 1, 1

    monkeypatch.setattr(combine, "_bigwig_to_bedgraph", fake_bedgraph)
    temp_dir = tmp_path / "tmp_bedgraphs"
    with pytest.raises(Exception):
        combine._combine_bigwig(
            [tmp_path / "chr1.bw"],
            tmp_path / "combined.bw",
            [("chr1", 10)],
            temp_dir,
        )
    assert temp_dir.exists()
    assert list(temp_dir.glob("*.bedGraph"))


def test_multicontig_manifest_can_defer_combination(tmp_path: Path) -> None:
    manifest = {
        "schema_version": 1,
        "command": "example",
        "combine_strategy": "native",
        "combined_name": "sample",
        "combined_dir": str(tmp_path / "combined"),
        "chrom_sizes": [["chr1", 10]],
        "per_contig": [],
    }
    path = tmp_path / combine.MANIFEST_NAME
    path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="no per-contig outputs"):
        combine._read_manifest(tmp_path)


def test_mnase_tree_combination_removes_worker_contig_from_filename(tmp_path):
    from nucleosuite.combine import combine_directory_trees

    chr1 = tmp_path / "per_contig" / "chr1" / "03_basic_tracks"
    chr2 = tmp_path / "per_contig" / "chr2" / "03_basic_tracks"
    chr1.mkdir(parents=True)
    chr2.mkdir(parents=True)
    (chr1 / "sample_chr1_fragments.bed").write_text("chr1\t1\t2\n")
    (chr2 / "sample_chr2_fragments.bed").write_text("chr2\t3\t4\n")

    output = tmp_path / "combined"
    combine_directory_trees(
        [chr1.parent, chr2.parent],
        output,
        chrom_sizes=[("chr1", 100), ("chr2", 100)],
        include_roots=["03_basic_tracks"],
        sample_name="sample",
    )

    combined = output / "03_basic_tracks" / "sample_fragments.bed"
    assert combined.read_text() == "chr1\t1\t2\nchr2\t3\t4\n"


def test_randomized_worker_names_keep_marker_terminal():
    from nucleosuite.cli.mnase_suite import _worker_sample_name

    assert (
        _worker_sample_name("sample_randomized_control", "chr1")
        == "sample_chr1_randomized_control"
    )


def test_multicontig_randomized_fragments_combine_qc_and_relocations(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from nucleosuite.cli.mnase_suite import _combine_randomized_fragments
    from nucleosuite import profile_plots

    roots = [tmp_path / "per_contig" / chrom for chrom in ("chr1", "chr2")]
    for index, (root, chrom) in enumerate(zip(roots, ("chr1", "chr2")), 1):
        setup = root / "00_setup"
        setup.mkdir(parents=True)
        prefix = setup / f"sample_{chrom}_randomized_control"
        with gzip.open(f"{prefix}.randomized.fragments.bed.gz", "wt") as handle:
            handle.write(f"{chrom}\t{index * 10}\t{index * 10 + 5}\n")
        Path(f"{prefix}.randomization_qc.tsv").write_text(
            "metric\tvalue\n"
            f"input\t{index}\n"
            f"matched\t{index}\n"
            "uniform\t0\n"
            "fallback\t0\n"
            "skipped\t0\n"
            f"unique_randomized_coordinates\t{index}\n"
            "duplicate_randomized_fragments\t0\n"
            f"maximum_randomized_multiplicity\t{index}\n"
            "collision_fraction\t0\n"
            "seed\t12345\n"
            "method\tdinucleotide\n"
        )
        Path(f"{prefix}.relocation_distances.tsv").write_text(
            "relocation_bp\tcount\n-10\t1\n10\t2\n"
        )

    monkeypatch.setattr(
        profile_plots,
        "plot_count_profile",
        lambda _input, output, **_kwargs: Path(output).write_bytes(b"png"),
    )
    combined_root = tmp_path / "combined"
    sample = "sample_randomized_control"
    bed = _combine_randomized_fragments(
        roots, combined_root, sample, ["chr1", "chr2"]
    )
    with gzip.open(bed, "rt") as handle:
        assert handle.read() == "chr1\t10\t15\nchr2\t20\t25\n"

    _fields, qc_rows = combine._read_tsv(
        combined_root / "00_setup" / f"{sample}.randomization_qc.tsv"
    )
    qc = {row["metric"]: row["value"] for row in qc_rows}
    assert qc["input"] == "3"
    assert qc["matched"] == "3"
    assert qc["maximum_randomized_multiplicity"] == "2"
    _fields, relocation_rows = combine._read_tsv(
        combined_root / "00_setup" / f"{sample}.relocation_distances.tsv"
    )
    assert {row["relocation_bp"]: row["count"] for row in relocation_rows} == {
        "-10": "2",
        "10": "4",
    }
    assert (
        combined_root / "00_setup" / f"{sample}.relocation_distances.png"
    ).read_bytes() == b"png"


def test_tree_combination_normalizes_bed_contigs(tmp_path):
    from nucleosuite.combine import combine_directory_trees

    one = tmp_path / "per_contig" / "20" / "01_combined_tracks"
    two = tmp_path / "per_contig" / "chr21" / "01_combined_tracks"
    one.mkdir(parents=True)
    two.mkdir(parents=True)
    (one / "sample_20_peaks.bed").write_text("20\t1\t2\n")
    (two / "sample_chr21_peaks.bed").write_text("chr21\t3\t4\n")
    output = tmp_path / "combined"
    combine_directory_trees(
        [one.parent, two.parent],
        output,
        chrom_sizes=[("chr20", 100), ("chr21", 100)],
        include_roots=["01_combined_tracks"],
        sample_name="sample",
    )
    combined = output / "01_combined_tracks" / "sample_peaks.bed"
    assert combined.read_text() == "chr20\t1\t2\nchr21\t3\t4\n"


def test_tree_combination_normalizes_terminal_randomized_marker(tmp_path):
    from nucleosuite.combine import combine_directory_trees

    one = tmp_path / "per_contig" / "chr1" / "01_combined_tracks"
    two = tmp_path / "per_contig" / "chr2" / "01_combined_tracks"
    one.mkdir(parents=True)
    two.mkdir(parents=True)
    (one / "sample_chr1_randomized_control_peaks.bed").write_text("chr1\t1\t2\n")
    (two / "sample_chr2_randomized_control_peaks.bed").write_text("chr2\t3\t4\n")
    output = tmp_path / "combined"
    combine_directory_trees(
        [one.parent, two.parent],
        output,
        chrom_sizes=[("chr1", 100), ("chr2", 100)],
        include_roots=["01_combined_tracks"],
        sample_name="sample_randomized_control",
    )
    combined = output / "01_combined_tracks" / "sample_randomized_control_peaks.bed"
    assert combined.read_text() == "chr1\t1\t2\nchr2\t3\t4\n"


def test_interval_combination_uses_bounded_fallback_sort_on_inversion(tmp_path):
    first = tmp_path / "chr1_a.bed"
    second = tmp_path / "chr2_b.bed"
    output = tmp_path / "combined.bed"
    first.write_text(
        "chr1\t20\t30\ttype1\n"
        "chr1\t10\t25\ttype2\n",
        encoding="utf-8",
    )
    second.write_text("chr2\t5\t15\ttype3\n", encoding="utf-8")

    combine._concatenate_intervals(
        [first, second],
        output,
        [("chr1", 100), ("chr2", 100)],
    )

    assert output.read_text(encoding="utf-8") == (
        "chr1\t10\t25\ttype2\n"
        "chr1\t20\t30\ttype1\n"
        "chr2\t5\t15\ttype3\n"
    )
    assert not list(tmp_path.glob(".nucleosuite_interval_sort_*"))


def test_interval_combination_preserves_already_sorted_order(tmp_path, monkeypatch):
    first = tmp_path / "chr1.bed"
    second = tmp_path / "chr2.bed"
    output = tmp_path / "combined.bed"
    first.write_text("chr1\t1\t2\nchr1\t3\t4\n", encoding="utf-8")
    second.write_text("chr2\t1\t2\n", encoding="utf-8")

    called = False

    def unexpected_sort(*_args, **_kwargs):
        nonlocal called
        called = True
        raise AssertionError("fallback sort should not run for ordered intervals")

    monkeypatch.setattr(combine, "_replace_with_sorted_interval_file", unexpected_sort)
    combine._concatenate_intervals(
        [first, second],
        output,
        [("chr1", 100), ("chr2", 100)],
    )
    assert called is False
    assert output.read_text(encoding="utf-8") == (
        "chr1\t1\t2\nchr1\t3\t4\nchr2\t1\t2\n"
    )


def test_bigwig_to_bedgraph_queries_bounded_chunks(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    class FakeHandle:
        def chroms(self):
            return {"chr1": 250_000}

        def intervals(self, chrom, start=None, end=None):
            assert chrom == "chr1"
            assert start is not None and end is not None
            assert end - start <= 100_000
            return [(start, min(start + 1, end), 1.0)]

        def close(self):
            pass

    class FakePyBigWig:
        @staticmethod
        def open(_path):
            return FakeHandle()

    monkeypatch.setattr(combine, "pyBigWig", FakePyBigWig())
    output = tmp_path / "track.bedGraph"
    intervals, chunks = combine._bigwig_to_bedgraph(
        tmp_path / "track.bw",
        output,
        ["chr1"],
        chrom_lengths={"chr1": 250_000},
        chunk_size=100_000,
    )
    assert intervals == 3
    assert chunks == 3
    assert output.read_text(encoding="utf-8").splitlines() == [
        "chr1\t0\t1\t1",
        "chr1\t100000\t100001\t1",
        "chr1\t200000\t200001\t1",
    ]


def test_tree_bedgraph_mode_consumes_validated_staged_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from nucleosuite.io.bedgraph import ValidatedBedGraphWriter

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    converter = bin_dir / "bedGraphToBigWig"
    _write_converter(converter, fail=False)
    monkeypatch.setenv("PATH", f"{bin_dir}:{Path('/usr/bin')}")
    monkeypatch.setattr(combine, "_verify_bigwig", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        combine,
        "_bigwig_to_bedgraph",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("staged mode must not reconvert BigWigs")
        ),
    )

    roots = []
    output = tmp_path / "combined"
    staging = output / "temporary_bedgraph_combine"
    for chrom, start, value in (("chr1", 0, 1.0), ("chr2", 5, 2.0)):
        root = tmp_path / "per_contig" / chrom
        track_dir = root / "01_combined_tracks"
        track_dir.mkdir(parents=True)
        bigwig = track_dir / f"sample_{chrom}_coverage.bw"
        bigwig.write_bytes(b"source-bigwig")
        staged = (
            staging
            / "per_contig"
            / chrom
            / "01_combined_tracks"
            / f"sample_{chrom}_coverage.bedGraph"
        )
        writer = ValidatedBedGraphWriter(
            staged,
            track="coverage",
            chrom_order=[chrom],
            source_bigwig=bigwig,
        )
        writer.add_interval(chrom, start, start + 2, value)
        writer.close()
        roots.append(root)

    result = combine.combine_directory_trees(
        roots,
        output,
        chrom_sizes=[("chr1", 100), ("chr2", 100)],
        include_roots=["01_combined_tracks"],
        sample_name="sample",
        bigwig_method="bedgraph",
    )

    combined = output / "01_combined_tracks" / "sample_coverage.bw"
    assert combined.read_text(encoding="utf-8") == (
        "chr1\t0\t2\t1\nchr2\t5\t7\t2\n"
    )
    assert result["bigwig_method"] == "bedgraph"
    assert not staging.exists()
    assert combine._combined_bigwig_marker_path(combined).is_file()


def test_direct_bigwig_combination_streams_bounded_chunks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = {
        str(tmp_path / "chr1.bw"): {"chr1": [(0, 2, 1.0), (150_000, 150_001, 2.0)]},
        str(tmp_path / "chr2.bw"): {"chr2": [(5, 7, 3.0)]},
    }
    written: dict[str, list[tuple[str, int, int, float]]] = {}
    headers: dict[str, list[tuple[str, int]]] = {}
    query_widths: list[int] = []

    for path in sources:
        Path(path).write_bytes(b"source")

    class ReadHandle:
        def __init__(self, path: str):
            self.path = path

        def chroms(self):
            return {chrom: 250_000 for chrom in sources[self.path]}

        def intervals(self, chrom, start=None, end=None):
            query_widths.append(int(end) - int(start))
            return [
                row for row in sources[self.path][chrom]
                if row[0] < end and row[1] > start
            ]

        def close(self):
            pass

    class WriteHandle:
        def __init__(self, path: str):
            self.path = path
            written[path] = []

        def addHeader(self, header):
            headers[self.path] = list(header)

        def addEntries(self, chroms, starts, ends=None, values=None, **_kwargs):
            if isinstance(chroms, str):
                chroms = [chroms] * len(starts)
            written[self.path].extend(
                (chrom, int(start), int(end), float(value))
                for chrom, start, end, value in zip(chroms, starts, ends, values)
            )

        def close(self):
            Path(self.path).write_bytes(b"combined")

    class VerifyHandle:
        def __init__(self, path: str):
            self.path = path

        def chroms(self):
            return dict(headers[self.path])

        def close(self):
            pass

    class FakePyBigWig:
        @staticmethod
        def open(path, mode=None):
            path = str(path)
            if mode == "w":
                return WriteHandle(path)
            if path in sources:
                return ReadHandle(path)
            return VerifyHandle(path)

    monkeypatch.setattr(combine, "pyBigWig", FakePyBigWig())
    output = tmp_path / "combined.bw"
    combine._combine_bigwig(
        [tmp_path / "chr1.bw", tmp_path / "chr2.bw"],
        output,
        [("chr1", 250_000), ("chr2", 250_000)],
        tmp_path / "unused",
        method="direct",
        chunk_size=100_000,
    )

    partial_key = str(output.with_name(output.name + ".partial"))
    assert headers[partial_key] == [("chr1", 250_000), ("chr2", 250_000)]
    assert written[partial_key] == [
        ("chr1", 0, 2, 1.0),
        ("chr1", 150_000, 150_001, 2.0),
        ("chr2", 5, 7, 3.0),
    ]
    assert max(query_widths) <= 100_000
    assert output.is_file()
    assert combine._combined_bigwig_marker_path(output).is_file()


def _write_track_report(path: Path, output_prefix: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "fragment_range\toutput_prefix\ttracks\tbasic_scope\tmax_duplicates\tmax_per_coordinate\n"
        f"120-180\t{output_prefix}\tpns,wps\tall\t1\t0\n",
        encoding="utf-8",
    )


def test_tree_combination_regenerates_one_combined_track_report(tmp_path: Path) -> None:
    from nucleosuite.combine import combine_directory_trees

    roots = [tmp_path / "per_contig" / name for name in ("chr1", "chr2")]
    for root in roots:
        prefix = root / "01_combined_tracks" / "pns" / f"sample_{root.name}_PNS"
        _write_track_report(
            root / "01_combined_tracks" / f"sample_{root.name}_completion_report.tsv",
            prefix,
        )

    output = tmp_path / "combined"
    result = combine_directory_trees(
        roots,
        output,
        chrom_sizes=[("chr1", 100), ("chr2", 100)],
        include_roots=["01_combined_tracks"],
        sample_name="sample",
        strict_complete=True,
    )

    assert result["warnings"] == []
    fields, rows = combine._read_tsv(
        output / "01_combined_tracks" / "sample_completion_report.tsv"
    )
    assert len(rows) == 1
    assert fields[1] == "output_prefix"
    assert rows[0]["output_prefix"] == str(
        output / "01_combined_tracks" / "pns" / "sample_PNS"
    )
    assert Path(result["combine_log"]).name == "combine_steps.log"


def test_tree_combination_rejects_missing_contig_contribution_in_strict_mode(
    tmp_path: Path,
) -> None:
    from nucleosuite.combine import combine_directory_trees

    roots = [tmp_path / "per_contig" / name for name in ("chr1", "chr2")]
    for root in roots:
        (root / "01_combined_tracks").mkdir(parents=True)
    (roots[0] / "01_combined_tracks" / "sample_chr1_peaks.bed").write_text(
        "chr1\t1\t2\n", encoding="utf-8"
    )

    with pytest.raises(ValueError, match="Refusing to create partial combined outputs"):
        combine_directory_trees(
            roots,
            tmp_path / "combined",
            chrom_sizes=[("chr1", 100), ("chr2", 100)],
            include_roots=["01_combined_tracks"],
            sample_name="sample",
            strict_complete=True,
        )


def test_tree_combination_rejects_unknown_analytical_tsv_in_strict_mode(
    tmp_path: Path,
) -> None:
    from nucleosuite.combine import combine_directory_trees

    roots = [tmp_path / "per_contig" / name for name in ("chr1", "chr2")]
    for root in roots:
        folder = root / "01_combined_tracks"
        folder.mkdir(parents=True)
        (folder / f"sample_{root.name}_analysis.tsv").write_text(
            "chromosome\tvalue\n" f"{root.name}\t1\n",
            encoding="utf-8",
        )

    with pytest.raises(ValueError, match="No explicit combination strategy"):
        combine_directory_trees(
            roots,
            tmp_path / "combined",
            chrom_sizes=[("chr1", 100), ("chr2", 100)],
            include_roots=["01_combined_tracks"],
            sample_name="sample",
            strict_complete=True,
        )


def test_directory_tree_combination_uses_requested_core_budget(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import threading
    import time

    roots = []
    for contig in ("chr1", "chr2"):
        root = tmp_path / "per_contig" / contig
        root.mkdir(parents=True)
        for name, start in (("a.bed", 1), ("b.bed", 5), ("c.bed", 9), ("d.bed", 13)):
            (root / name).write_text(f"{contig}\t{start}\t{start + 1}\n", encoding="utf-8")
        roots.append(root)

    original = combine._concatenate_intervals
    lock = threading.Lock()
    active = 0
    maximum_active = 0

    def measured(inputs, output, chrom_sizes):
        nonlocal active, maximum_active
        with lock:
            active += 1
            maximum_active = max(maximum_active, active)
        try:
            time.sleep(0.05)
            return original(inputs, output, chrom_sizes)
        finally:
            with lock:
                active -= 1

    monkeypatch.setattr(combine, "_concatenate_intervals", measured)
    result = combine.combine_directory_trees(
        roots,
        tmp_path / "combined",
        chrom_sizes=[("chr1", 100), ("chr2", 100)],
        cores=2,
    )

    assert result["cores"] == 2
    assert maximum_active == 2
    for name in ("a.bed", "b.bed", "c.bed", "d.bed"):
        assert (tmp_path / "combined" / name).exists()


def test_tracks_and_combine_parsers_expose_core_budget() -> None:
    from nucleosuite.cli.main import build_parser

    parser = build_parser()
    tracks = parser.parse_args(
        ["tracks", "--bam", "sample.bam", "--fragment-range", "147=dyad", "--cores", "4"]
    )
    assert tracks.cores == 4

    combined = combine.build_parser().parse_args(["--input-dir", "run", "--cores", "3"])
    assert combined.cores == 3
