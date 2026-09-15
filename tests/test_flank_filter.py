from __future__ import annotations

from pathlib import Path

from nucleosuite.flank_filter import (
    build_flank_index,
    filter_flanked_regions,
    main,
    read_bed_records,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def test_flank_filter_keeps_only_regions_with_both_sides(tmp_path: Path):
    regions = _write(
        tmp_path / "regions.bed",
        "chr1\t90\t110\tleft_only\n"
        "chr1\t190\t210\tflanked\n"
        "chr1\t290\t310\tright_only\n"
        "chr2\t190\t210\twrong_chrom\n",
    )
    flanks = _write(
        tmp_path / "flanks.bed",
        "chr1\t145\t155\n"
        "chr1\t245\t255\n"
        "chr2\t245\t255\n",
    )

    region_rows = read_bed_records(regions)
    retained = filter_flanked_regions(
        region_rows,
        build_flank_index(read_bed_records(flanks)),
    )
    assert [record.raw_line.split("\t")[3] for record in retained] == ["flanked"]


def test_default_position_prefers_bed_column_7_and_preserves_rows(tmp_path: Path):
    regions = _write(
        tmp_path / "regions.bed",
        # Midpoint is 500 but column 7 is 100. The flanks surround 100 only.
        "chr1\t0\t1000\tnuc1\t12.5\t.\t100\t101\textra\n",
    )
    flanks = _write(
        tmp_path / "breakpoints.bed",
        "chr1\t40\t60\tbp1\t-5\t.\t50\t51\n"
        "chr1\t140\t160\tbp2\t-4\t.\t150\t151\n",
    )
    output = tmp_path / "filtered.bed"

    assert main(["--regions", str(regions), "--flanks", str(flanks), "--out", str(output)]) == 0
    assert output.read_text() == regions.read_text()


def test_bed3_falls_back_to_interval_midpoint(tmp_path: Path):
    regions = _write(tmp_path / "regions.bed", "chr1\t90\t110\n")
    flanks = _write(tmp_path / "flanks.bed", "chr1\t40\t60\nchr1\t140\t160\n")
    retained = filter_flanked_regions(
        read_bed_records(regions),
        build_flank_index(read_bed_records(flanks)),
    )
    assert len(retained) == 1
    assert retained[0].position == 100


def test_equal_position_is_not_counted_as_either_side(tmp_path: Path):
    regions = _write(tmp_path / "regions.bed", "chr1\t90\t110\n")
    flanks = _write(tmp_path / "flanks.bed", "chr1\t90\t110\nchr1\t140\t160\n")
    retained = filter_flanked_regions(
        read_bed_records(regions),
        build_flank_index(read_bed_records(flanks)),
    )
    assert retained == []


def test_max_flank_distance_requires_both_nearest_flanks_within_limit(tmp_path: Path):
    regions = _write(tmp_path / "regions.bed", "chr1\t95\t105\n")
    flanks = _write(tmp_path / "flanks.bed", "chr1\t45\t55\nchr1\t195\t205\n")
    region_rows = read_bed_records(regions)
    index = build_flank_index(read_bed_records(flanks))

    assert len(filter_flanked_regions(region_rows, index)) == 1
    assert filter_flanked_regions(region_rows, index, max_flank_distance=75) == []
    assert len(filter_flanked_regions(region_rows, index, max_flank_distance=100)) == 1


def test_explicit_position_columns_override_auto_bed7(tmp_path: Path):
    regions = _write(
        tmp_path / "regions.bed",
        "chr1\t0\t1000\tnuc1\t0\t.\t900\t901\t100\n",
    )
    flanks = _write(
        tmp_path / "flanks.bed",
        "chr1\t0\t10\tbp1\t0\t.\t800\t801\t50\n"
        "chr1\t0\t10\tbp2\t0\t.\t850\t851\t150\n",
    )
    retained = filter_flanked_regions(
        read_bed_records(regions, position_column=9),
        build_flank_index(read_bed_records(flanks, position_column=9)),
    )
    assert len(retained) == 1
    assert retained[0].position == 100


def test_consecutive_regions_between_flanks_are_both_rejected(tmp_path: Path):
    regions = _write(
        tmp_path / "regions.bed",
        "chr1\t95\t105\tnuc_a\n"
        "chr1\t115\t125\tnuc_b\n",
    )
    flanks = _write(
        tmp_path / "flanks.bed",
        "chr1\t45\t55\tbrk_left\n"
        "chr1\t165\t175\tbrk_right\n",
    )

    retained = filter_flanked_regions(
        read_bed_records(regions),
        build_flank_index(read_bed_records(flanks)),
    )
    assert retained == []


def test_consecutive_flanks_do_not_disqualify_a_clean_region(tmp_path: Path):
    regions = _write(tmp_path / "regions.bed", "chr1\t195\t205\tnuc\n")
    flanks = _write(
        tmp_path / "flanks.bed",
        "chr1\t45\t55\tbrk_1\n"
        "chr1\t145\t155\tbrk_2\n"
        "chr1\t245\t255\tbrk_3\n",
    )

    retained = filter_flanked_regions(
        read_bed_records(regions),
        build_flank_index(read_bed_records(flanks)),
    )
    assert [record.raw_line.split("\t")[3] for record in retained] == ["nuc"]


def test_region_outside_previous_flank_pair_does_not_disqualify_clean_region(tmp_path: Path):
    regions = _write(
        tmp_path / "regions.bed",
        "chr1\t45\t55\tleft_nuc\n"
        "chr1\t195\t205\ttarget_nuc\n",
    )
    flanks = _write(
        tmp_path / "flanks.bed",
        "chr1\t145\t155\tbrk_left\n"
        "chr1\t245\t255\tbrk_right\n",
    )

    retained = filter_flanked_regions(
        read_bed_records(regions),
        build_flank_index(read_bed_records(flanks)),
    )
    assert [record.raw_line.split("\t")[3] for record in retained] == ["target_nuc"]


def test_duplicate_region_positions_are_rejected(tmp_path: Path):
    regions = _write(
        tmp_path / "regions.bed",
        "chr1\t95\t105\tnuc_a\n"
        "chr1\t95\t105\tnuc_b\n",
    )
    flanks = _write(
        tmp_path / "flanks.bed",
        "chr1\t45\t55\tbrk_left\n"
        "chr1\t145\t155\tbrk_right\n",
    )

    retained = filter_flanked_regions(
        read_bed_records(regions),
        build_flank_index(read_bed_records(flanks)),
    )
    assert retained == []
