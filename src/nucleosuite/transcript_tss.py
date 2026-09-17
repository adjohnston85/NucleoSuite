"""Extract gene-matched transcript starts from Ensembl GRCh37 release 87."""

from __future__ import annotations

import csv
import gzip
import os
import re
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.request import urlopen

from nucleosuite.io import open_text

ENSEMBL_GTF_URL = (
    "https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/"
    "Homo_sapiens.GRCh37.87.gtf.gz"
)
TSS_FILENAME = "hg19_ensembl87_transcript_tss.tsv.gz"
HEADER = ("gene_id", "transcript_id", "chrom", "tss_start", "tss_end", "strand")
_ATTRIBUTE = re.compile(r'(gene_id|transcript_id)\s+"([^"]+)"')


@dataclass(frozen=True)
class TranscriptTSS:
    gene_id: str
    transcript_id: str
    chrom: str
    start: int
    end: int
    strand: str


def _base_id(value: str) -> str:
    return value.split(".", 1)[0]


def _match_chrom(chrom: str, expected: str) -> bool:
    return chrom.removeprefix("chr") == expected.removeprefix("chr")


def extract_transcript_tss(gtf: str | Path, genes, output: str | Path) -> tuple[int, int]:
    """Write BED-coordinate transcript TSSs for matching Ensembl gene IDs.

    The input GTF is one-based inclusive; output TSS intervals are zero-based,
    half-open and one base long. Only genuine transcript records are used.
    """
    by_id = {_base_id(gene.gene_id): gene for gene in genes}
    found: dict[tuple[str, str], TranscriptTSS] = {}
    with open_text(gtf) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            attrs = dict(_ATTRIBUTE.findall(parts[8]))
            gene_id = _base_id(attrs.get("gene_id", ""))
            transcript_id = _base_id(attrs.get("transcript_id", ""))
            gene = by_id.get(gene_id)
            if not gene or not transcript_id or parts[6] not in ("+", "-"):
                continue
            if not _match_chrom(parts[0], gene.chrom) or len(gene.fields) < 6 or parts[6] != gene.fields[5]:
                continue
            start, end = int(parts[3]), int(parts[4])
            tss = start - 1 if parts[6] == "+" else end - 1
            if not (gene.start <= tss < gene.end):
                continue
            record = TranscriptTSS(gene.gene_id, transcript_id, gene.chrom, tss, tss + 1, parts[6])
            key = (gene.gene_id, transcript_id)
            previous = found.get(key)
            if previous is not None and previous != record:
                raise ValueError(f"Conflicting coordinates for transcript {transcript_id}")
            found[key] = record
    if not found:
        raise ValueError("No transcript TSSs from the GTF match the gene BED (check assembly and Ensembl release)")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", encoding="utf-8", newline="") if output.suffix == ".gz" else output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        for key in sorted(found, key=lambda key: (found[key].chrom, found[key].start, key)):
            record = found[key]
            writer.writerow((record.gene_id, record.transcript_id, record.chrom, record.start, record.end, record.strand))
    return len(found), len({item.gene_id for item in found.values()})


def read_transcript_tss(path: str | Path, genes) -> dict[str, tuple[TranscriptTSS, ...]]:
    """Validate all TSS rows, including Ensembl ID and strand matching."""
    gene_index = {_base_id(gene.gene_id): gene for gene in genes}
    by_gene: dict[str, dict[str, TranscriptTSS]] = defaultdict(dict)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(HEADER).issubset(reader.fieldnames):
            raise ValueError(f"Transcript annotation requires TSV columns: {', '.join(HEADER)}")
        for line_num, row in enumerate(reader, 2):
            gene = gene_index.get(_base_id(row["gene_id"]))
            if gene is None:
                continue
            try:
                start, end = int(row["tss_start"]), int(row["tss_end"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_num}: invalid transcript TSS coordinates") from exc
            if (end != start + 1 or not gene.start <= start < gene.end
                    or not _match_chrom(row["chrom"], gene.chrom)
                    or len(gene.fields) < 6 or row["strand"] != gene.fields[5]
                    or not row["transcript_id"]):
                raise ValueError(f"{path}:{line_num}: TSS does not match the gene interval, strand or assembly")
            item = TranscriptTSS(gene.gene_id, row["transcript_id"], gene.chrom, start, end, row["strand"])
            prior = by_gene[gene.gene_id].get(item.transcript_id)
            if prior is not None and prior != item:
                raise ValueError(f"{path}:{line_num}: transcript ID has conflicting TSS coordinates")
            by_gene[gene.gene_id][item.transcript_id] = item
    if not by_gene:
        raise ValueError(f"No transcript TSSs match the supplied genes in {path}")
    return {gene_id: tuple(sorted(entries.values(), key=lambda t: (t.start, t.transcript_id)))
            for gene_id, entries in by_gene.items()}


def ensembl_cache_path() -> Path:
    root = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return root / "nucleosuite" / TSS_FILENAME


def ensure_ensembl_tss(genes) -> Path:
    """Create a user-writable cached transcript annotation from Ensembl once."""
    dest = ensembl_cache_path()
    if dest.is_file() and dest.stat().st_size > 64:
        existing = read_transcript_tss(dest, genes)
        if len(existing) >= max(1, int(len(genes) * 0.8)):
            return dest
        # Recreate an annotation that was cached for a smaller gene subset.
    dest.parent.mkdir(parents=True, exist_ok=True)
    source = None
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".gtf.gz", dir=dest.parent, delete=False) as handle:
            source = Path(handle.name)
            try:
                with urlopen(ENSEMBL_GTF_URL, timeout=90) as response:
                    while chunk := response.read(1024 * 1024):
                        handle.write(chunk)
            except OSError as exc:
                raise RuntimeError(
                    "Could not obtain the Ensembl GRCh37 release-87 transcript annotation. "
                    "Download Homo_sapiens.GRCh37.87.gtf.gz from " + ENSEMBL_GTF_URL +
                    " and extract transcript TSSs with examples/build_ensembl87_transcript_tss.py; "
                    "then pass --transcript-tss-tsv FILE."
                ) from exc
        with tempfile.NamedTemporaryFile(suffix=".tsv.gz", dir=dest.parent, delete=False) as handle:
            temporary = Path(handle.name)
        _, matching = extract_transcript_tss(source, genes, temporary)
        if matching < max(1, int(len(genes) * 0.8)):
            raise ValueError(f"Ensembl annotation matches only {matching}/{len(genes)} genes; check annotation provenance")
        temporary.replace(dest)
        return dest
    finally:
        if source is not None:
            source.unlink(missing_ok=True)
        if temporary is not None:
            temporary.unlink(missing_ok=True)
