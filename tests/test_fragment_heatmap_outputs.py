from __future__ import annotations

from pathlib import Path

from nucleosuite.fragment_heatmap import main


def _profiles(tmp_path: Path) -> tuple[Path, Path]:
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    a.write_text("fragment_length\tcount\n100\t10\n101\t20\n102\t10\n", encoding="utf-8")
    b.write_text("fragment_length\tcount\n100\t20\n101\t10\n102\t20\n", encoding="utf-8")
    return a, b


def test_fragment_heatmap_keeps_compact_matrix_but_omits_workbook_by_default(tmp_path: Path) -> None:
    a, b = _profiles(tmp_path)
    prefix = tmp_path / "heatmap"
    assert main([
        "--input", f"A={a}", "--input", f"B={b}", "--out-prefix", str(prefix),
        "--min-frag", "100", "--max-frag", "102", "--no-cluster",
    ]) == 0
    assert list(tmp_path.glob("heatmap*_normalised_matrix.tsv"))
    assert list(tmp_path.glob("heatmap*_heatmap.png"))
    assert not list(tmp_path.glob("heatmap*.xlsx"))


def test_fragment_heatmap_workbook_is_opt_in_detail_output(tmp_path: Path) -> None:
    a, b = _profiles(tmp_path)
    prefix = tmp_path / "heatmap_detail"
    assert main([
        "--input", f"A={a}", "--input", f"B={b}", "--out-prefix", str(prefix),
        "--min-frag", "100", "--max-frag", "102", "--no-cluster", "--write-detail-tables",
    ]) == 0
    assert list(tmp_path.glob("heatmap_detail*_normalised_matrix.tsv"))
    assert list(tmp_path.glob("heatmap_detail*.xlsx"))
