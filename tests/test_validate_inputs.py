from __future__ import annotations

import gzip


def test_validate_interval_reads_gzip_to_eof_and_reports_corruption(tmp_path):
    from nucleosuite.validate_inputs import _validate_interval_file

    good = tmp_path / "fragments.bed.gz"
    with gzip.open(good, "wt") as handle:
        handle.write("chr1\t1\t5\nchr1\t8\t12\n")
    row = _validate_interval_file(
        good, kind="fragments", max_records=None, require_sorted=True
    )
    assert row.status == "PASS"
    assert row.records == 2

    corrupt = tmp_path / "corrupt.bed.gz"
    corrupt.write_bytes(good.read_bytes()[:-6])
    row = _validate_interval_file(
        corrupt, kind="fragments", max_records=None, require_sorted=False
    )
    assert row.status == "FAIL"


def test_reference_compatibility_resolves_conservative_aliases():
    from nucleosuite.validate_inputs import _reference_compatibility

    row = _reference_compatibility(
        [("bam", {"1": 100, "MT": 20}), ("fasta", {"chr1": 100, "chrM": 20})]
    )
    assert row is not None and row.status == "PASS"

    row = _reference_compatibility(
        [("bam", {"1": 100}), ("fasta", {"chr1": 101})]
    )
    assert row is not None and row.status == "FAIL"


def test_validation_report_is_atomic_tsv(tmp_path):
    from nucleosuite.validate_inputs import ValidationRow, _write_report

    report = tmp_path / "validation.tsv"
    _write_report(report, [ValidationRow("bed", "input.bed", "PASS", 1, "ok")])
    assert report.read_text().splitlines() == [
        "kind\tpath\tstatus\trecords\tdetail",
        "bed\tinput.bed\tPASS\t1\tok",
    ]
    assert not list(tmp_path.glob("*.partial"))
