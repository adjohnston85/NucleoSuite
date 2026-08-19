"""Tests for gene-level peak-spacing and FFT/expression analyses."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import sys
from types import SimpleNamespace

from nucleosuite import gene_expression



class FakeBigWig:
    def __init__(self, values: np.ndarray):
        self._values = np.asarray(values, dtype=float)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False

    def chroms(self):
        return {"chr1": len(self._values)}

    def values(self, chrom, start, end, numpy=False):
        assert chrom == "chr1"
        result = self._values[start:end].copy()
        return result if numpy else result.tolist()


def test_periodogram_reports_dominant_period():
    positions = np.arange(1000, dtype=float)
    values = np.sin(2.0 * np.pi * positions / 50.0)
    periods = np.arange(40, 81)
    intensities = gene_expression.periodogram_intensities(
        values,
        periods,
        trim_fraction=0.1,
        pad_fraction=0.3,
        taper_fraction=0.3,
        recursive=False,
    )
    dominant = int(periods[int(np.argmax(intensities))])
    assert 48 <= dominant <= 52


def test_all_valid_fft_mask_matches_unmasked_numerics():
    positions = np.arange(600, dtype=float)
    values = np.sin(2.0 * np.pi * positions / 55.0)
    periods = np.arange(40, 81)
    kwargs = dict(
        trim_fraction=0.1,
        pad_fraction=0.3,
        taper_fraction=0.3,
        recursive=True,
    )
    unmasked = gene_expression.periodogram_intensities(values, periods, **kwargs)
    masked = gene_expression.periodogram_intensities(
        values, periods, valid_mask=np.ones(values.size, dtype=bool), **kwargs
    )
    assert np.array_equal(masked, unmasked)


def test_gene_blacklist_filter_uses_only_strand_aware_tss_anchors(tmp_path):
    from nucleosuite.core.blacklist import load_blacklist_unbounded

    genes = [
        gene_expression.GeneRecord("chr1", 10, 100, "plus_body", "plus_body", "+"),
        gene_expression.GeneRecord("chr1", 120, 200, "minus_tss", "minus_tss", "-"),
    ]
    blacklist_path = tmp_path / "blacklist.bed"
    blacklist_path.write_text("chr1\t50\t60\nchr1\t199\t200\n")
    retained, excluded = gene_expression.filter_blacklisted_gene_anchors(
        genes, load_blacklist_unbounded(blacklist_path)
    )
    assert [gene.gene_id for gene in retained] == ["plus_body"]
    assert excluded == 1


def test_gene_expression_all_outputs(tmp_path, monkeypatch):
    genes = tmp_path / "genes.bed"
    genes.write_text(
        "chr1\t100\t1100\tENSG00000000001\tGENE1\t+\n"
        "chr1\t1300\t2300\tENSG00000000002\tGENE2\t+\n"
        "chr1\t2500\t3500\tENSG00000000003\tGENE3\t+\n"
    )

    expression = tmp_path / "expression.tsv"
    with expression.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Gene", "Gene name", "Cell line", "TPM", "pTPM", "nTPM"])
        for gene, name, base in (
            ("ENSG00000000001", "GENE1", 1.0),
            ("ENSG00000000002", "GENE2", 5.0),
            ("ENSG00000000003", "GENE3", 20.0),
        ):
            for profile, multiplier in (("NB-4", 1.0), ("CellB", 2.0), ("CellC", 0.5)):
                value = base * multiplier
                writer.writerow([gene, name, profile, value, value, value])

    peaks = tmp_path / "peaks.bed"
    rows = []
    for gene_index, start in enumerate((100, 1300, 2500)):
        spacing = (200, 175, 150)[gene_index]
        for peak_index in range(5):
            centre = start + 50 + peak_index * spacing
            rows.append(
                ["chr1", centre - 20, centre + 21, f"p{gene_index}_{peak_index}", 10, ".", centre, 10]
            )
    peaks.write_text("".join("\t".join(map(str, row)) + "\n" for row in rows))

    x = np.arange(5000, dtype=float)
    signal = np.zeros_like(x)
    for start, period in ((100, 60), (1300, 55), (2500, 50)):
        signal[start : start + 1000] = np.sin(2 * np.pi * np.arange(1000) / period)
    bigwig = tmp_path / "pns.bw"
    bigwig.touch()
    monkeypatch.setitem(sys.modules, "pyBigWig", SimpleNamespace(open=lambda _: FakeBigWig(signal)))

    prefix = tmp_path / "result" / "sample"
    parser = gene_expression.build_parser()
    args = parser.parse_args(
        [
            "--expression", str(expression),
            "--genes-bed", str(genes),
            "--peaks", f"sample={peaks}",
            "--signal", f"sample={bigwig}",
            "--analysis", "all",
            "--output-prefix", str(prefix),
            "--gene-flank", "100",
            "--high-confidence-peaks", "3",
            "--min-correlation-genes", "3",
            "--fft-window", "500",
            "--fft-period-min", "40",
            "--fft-period-max", "80",
            "--fft-ranking-periods", "50,55,60",
            "--no-fft-recursive-filter",
            "--top-profiles", "3",
            "--focus-profile", "NB-4",
        ]
    )
    assert gene_expression.run(args) == 0
    prefix = Path(
        f"{prefix}_analysisall_flank100_fftwin500"
    )

    expected = [
        f"{prefix}_gene_peak_spacing.tsv",
        f"{prefix}_spacing_expression_correlations.tsv",
        f"{prefix}_spacing_expression_scatter.png",
        f"{prefix}_per_gene_fft.tsv.gz",
        f"{prefix}_fft_expression_correlations.tsv",
        f"{prefix}_expression_profile_rankings.tsv",
        f"{prefix}_metadata.tsv",
    ]
    for path in expected:
        assert Path(path).is_file()
        assert Path(path).stat().st_size > 0

    spacing_rows = list(csv.DictReader(Path(f"{prefix}_gene_peak_spacing.tsv").open(), delimiter="\t"))
    assert [float(row["body_median_spacing_bp"]) for row in spacing_rows] == [200.0, 175.0, 150.0]


def test_combined_profiles_use_profile_type_for_cell_line_metadata(tmp_path, monkeypatch):
    genes = tmp_path / "genes.bed"
    genes.write_text(
        "chr1\t100\t1100\tENSG00000000001\tGENE1\t+\n"
        "chr1\t1300\t2300\tENSG00000000002\tGENE2\t+\n"
        "chr1\t2500\t3500\tENSG00000000003\tGENE3\t+\n"
    )

    expression = tmp_path / "combined_expression.tsv"
    with expression.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["Gene", "Gene name", "Profile", "Profile type", "Cell line", "Tissue", "TPM", "pTPM", "nTPM"]
        )
        for gene_index, (gene, name) in enumerate(
            (
                ("ENSG00000000001", "GENE1"),
                ("ENSG00000000002", "GENE2"),
                ("ENSG00000000003", "GENE3"),
            ),
            1,
        ):
            writer.writerow([gene, name, "CellA", "cell_line", "CellA", "", gene_index, gene_index, gene_index])
            writer.writerow([gene, name, "bone marrow", "tissue", "", "bone marrow", gene_index, gene_index, gene_index])

    metadata = tmp_path / "cell_lines.tsv"
    metadata.write_text(
        "Cell line\tCancer cell line\tPrimary disease\tDisease subtype\tPrimary/Metastasis\tSample collection site\tCellosaurus ID\n"
        "CellA\tYes\tLeukaemia\tExample subtype\tPrimary\tBlood\tCVCL_0001\n"
    )

    signal_path = tmp_path / "signal.bw"
    signal_path.touch()
    x = np.arange(5000, dtype=float)
    signal = np.sin(2 * np.pi * x / 50.0)
    monkeypatch.setitem(sys.modules, "pyBigWig", SimpleNamespace(open=lambda _: FakeBigWig(signal)))

    prefix = tmp_path / "result" / "sample"
    args = gene_expression.build_parser().parse_args(
        [
            "--expression", str(expression),
            "--genes-bed", str(genes),
            "--signal", f"sample={signal_path}",
            "--analysis", "fft",
            "--output-prefix", str(prefix),
            "--expression-profile-column", "Profile",
            "--expression-value-column", "nTPM",
            "--profile-metadata", str(metadata),
            "--min-nonzero-profiles-per-gene", "1",
            "--min-correlation-genes", "3",
            "--fft-window", "500",
            "--fft-period-min", "40",
            "--fft-period-max", "80",
            "--fft-ranking-periods", "50",
            "--no-fft-recursive-filter",
        ]
    )
    assert gene_expression.run(args) == 0
    prefix = Path(
        f"{prefix}_analysisfft_flank10000_fftwin500"
    )

    rows = list(
        csv.DictReader(
            Path(f"{prefix}_expression_profile_rankings.tsv").open(),
            delimiter="\t",
        )
    )
    by_profile = {row["profile"]: row for row in rows}
    assert by_profile["CellA"]["profile_type"] == "cell_line"
    assert by_profile["CellA"]["metadata_matched"] == "true"
    assert by_profile["CellA"]["metadata_status"] == "matched"
    assert by_profile["CellA"]["cellosaurus_id"] == "CVCL_0001"
    assert by_profile["bone marrow"]["profile_type"] == "tissue"
    assert by_profile["bone marrow"]["metadata_matched"] == "false"
    assert by_profile["bone marrow"]["metadata_status"] == "not_applicable"
