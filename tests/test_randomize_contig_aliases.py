from __future__ import annotations

import gzip
from types import SimpleNamespace

import pytest


@pytest.mark.parametrize(
    ("source_contig", "fasta_contig"),
    [("21", "chr21"), ("chr21", "21")],
)
def test_randomize_fragments_resolves_bam_fasta_contig_aliases(
    tmp_path, monkeypatch, source_contig, fasta_contig
):
    import nucleosuite.randomize_fragments_command as command

    sequence = "A" * 500
    fetched_names: list[str] = []

    class FakeFasta:
        references = (fasta_contig,)

        def get_reference_length(self, name):
            assert name == fasta_contig
            return len(sequence)

        def fetch(self, name, start, end):
            fetched_names.append(name)
            assert name == fasta_contig
            return sequence[start:end]

        def close(self):
            return None

    class FakeSource:
        references = [source_contig]
        lengths = [len(sequence)]

        def fetch(self, contig, start, end, **_kwargs):
            assert contig == source_contig
            return [(10, 155)]

        def close(self):
            return None

    monkeypatch.setattr(
        command,
        "pysam",
        SimpleNamespace(FastaFile=lambda _path: FakeFasta()),
    )
    monkeypatch.setattr(command, "open_fragment_source", lambda **_kwargs: FakeSource())
    monkeypatch.setattr(command, "plot_count_profile", lambda *_args, **_kwargs: None)

    prefix = tmp_path / "randomized"
    args = command.build_parser().parse_args(
        [
            "--bam",
            "input.bam",
            "--fasta",
            "reference.fa",
            "--output-prefix",
            str(prefix),
            "--contigs",
            source_contig,
            "--frag-lower",
            "145",
            "--frag-upper",
            "145",
            "--search-window",
            "500",
        ]
    )

    assert command._run_serial(args) == 0
    output = tmp_path / "randomized.randomized.fragments.bed.gz"
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        fields = handle.readline().rstrip("\n").split("\t")
    assert fields[0] == source_contig
    assert int(fields[2]) - int(fields[1]) == 145
    assert fetched_names and set(fetched_names) == {fasta_contig}


def test_randomize_fragments_rejects_alias_length_mismatch(tmp_path, monkeypatch):
    import nucleosuite.randomize_fragments_command as command

    class FakeFasta:
        references = ("chr21",)

        def get_reference_length(self, _name):
            return 499

        def close(self):
            return None

    class FakeSource:
        references = ["21"]
        lengths = [500]

        def close(self):
            return None

    monkeypatch.setattr(
        command,
        "pysam",
        SimpleNamespace(FastaFile=lambda _path: FakeFasta()),
    )
    monkeypatch.setattr(command, "open_fragment_source", lambda **_kwargs: FakeSource())

    args = command.build_parser().parse_args(
        [
            "--bam",
            "input.bam",
            "--fasta",
            "reference.fa",
            "--output-prefix",
            str(tmp_path / "randomized"),
            "--contigs",
            "21",
            "--search-window",
            "500",
        ]
    )
    with pytest.raises(SystemExit) as exc:
        command._run_serial(args)
    assert exc.value.code == 2
