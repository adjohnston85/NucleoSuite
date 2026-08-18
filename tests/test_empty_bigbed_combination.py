from __future__ import annotations

from pathlib import Path


def test_empty_bigbed_markers_count_as_complete_contributions(tmp_path):
    from nucleosuite.combine import combine_directory_trees

    roots = [tmp_path / "chr1", tmp_path / "chr2"]
    relative = Path("01_combined_tracks/sequence/unclassified.bb.empty")
    for root in roots:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("status\tempty\nrecord_count\t0\n", encoding="utf-8")

    output = tmp_path / "combined"
    result = combine_directory_trees(
        roots,
        output,
        chrom_sizes=[("chr1", 100), ("chr2", 100)],
        strict_complete=True,
    )
    assert result["warnings"] == []
    assert (output / "01_combined_tracks/sequence/unclassified.bed").exists()
    marker = output / "01_combined_tracks/sequence/unclassified.bb.empty"
    assert marker.exists()
    assert not (output / "01_combined_tracks/sequence/unclassified.bb").exists()


def test_mixed_empty_and_nonempty_bigbed_is_not_incomplete(tmp_path, monkeypatch):
    import nucleosuite.combine as combine

    root1 = tmp_path / "chr1"
    root2 = tmp_path / "chr2"
    relative = Path("01_combined_tracks/sequence/type1.bb")
    actual = root1 / relative
    actual.parent.mkdir(parents=True, exist_ok=True)
    actual.write_bytes(b"placeholder")
    marker = root2 / Path(str(relative) + ".empty")
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text("status\tempty\nrecord_count\t0\n", encoding="utf-8")

    def fake_combine(inputs, output_bed, chrom_sizes=None):
        assert inputs == [actual]
        Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
        Path(output_bed).write_text("chr1\t1\t2\n", encoding="utf-8")

    monkeypatch.setattr(combine, "_combine_bigbed", fake_combine)
    monkeypatch.setattr(combine.shutil, "which", lambda _name: None)
    result = combine.combine_directory_trees(
        [root1, root2],
        tmp_path / "combined",
        chrom_sizes=[("chr1", 100), ("chr2", 100)],
        strict_complete=True,
    )
    assert not any("missing" in warning.lower() for warning in result["warnings"])
