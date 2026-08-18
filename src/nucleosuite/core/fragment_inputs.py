"""Unified fragment-coordinate input for BAM, BED, BED.gz and bigBed files.

BAM headers define a canonical output namespace.  Materialised fragment files
are stored in a temporary on-disk SQLite index so large BED/BED.gz collections
can be queried chunk by chunk without being loaded into memory.
"""

from __future__ import annotations

import os
import random
import sqlite3
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterator, Sequence

try:
    import pysam
except ImportError:  # allow help/tests without optional runtime dependency
    pysam = None

try:
    import pyBigWig
except ImportError:  # allow help/tests without optional runtime dependency
    pyBigWig = None

from nucleosuite.core.bam_headers import merge_bam_reference_headers_with_aliases
from nucleosuite.core.blacklist import BlacklistIndex, load_blacklist
from nucleosuite.core.chrom_sizes import chromosome_size_dict
from nucleosuite.core.fragments import collect_fragments, require_bam_indexes
from nucleosuite.core.regions import canonical_contig_key, resolve_contig_name
from nucleosuite.io import open_text, strip_known_suffix

Fragment = tuple[int, int]


def _read_chrom_sizes(path: str | os.PathLike[str] | None) -> dict[str, int]:
    return chromosome_size_dict(path)


def _choose_canonical_name(names: Sequence[str]) -> str:
    for name in names:
        if name.lower().startswith("chr"):
            return name
    return names[0]


class FragmentSource:
    """Abstract fetchable fragment source."""

    references: list[str]
    lengths: list[int]
    input_paths: list[str]
    kind: str

    def fetch(
        self,
        contig: str,
        start: int,
        end: int,
        *,
        max_per_coordinate: int,
        subsample: float | None,
        dedup_scope: str,
    ) -> list[Fragment]:
        raise NotImplementedError

    def close(self) -> None:
        return None

    @property
    def label(self) -> str:
        stems = [strip_known_suffix(path) for path in self.input_paths]
        return stems[0] if len(stems) == 1 else "_".join(stems)




class BlacklistFragmentSource(FragmentSource):
    """Filter complete fragments that overlap an optional BED blacklist."""

    def __init__(self, source: FragmentSource, blacklist: BlacklistIndex):
        self.source = source
        self.blacklist = blacklist
        self.references = source.references
        self.lengths = source.lengths
        self.input_paths = source.input_paths
        self.kind = source.kind
        self.fragments_examined = 0
        self.fragments_excluded = 0

    def fetch(
        self, contig: str, start: int, end: int, *,
        max_per_coordinate: int, subsample: float | None, dedup_scope: str,
    ) -> list[Fragment]:
        fragments = self.source.fetch(
            contig, start, end, max_per_coordinate=max_per_coordinate,
            subsample=subsample, dedup_scope=dedup_scope,
        )
        canonical = resolve_contig_name(
            contig, self.references, source_label="fragment source"
        )
        self.fragments_examined += len(fragments)
        retained = [
            fragment for fragment in fragments
            if not self.blacklist.overlaps(canonical, fragment[0], fragment[1])
        ]
        self.fragments_excluded += len(fragments) - len(retained)
        return retained

    def close(self) -> None:
        self.source.close()

class BamFragmentSource(FragmentSource):
    """Fetch paired-end fragments using BAM-derived canonical contig names."""

    kind = "bam"

    def __init__(self, paths: Sequence[str], *, chrom_sizes: str | None = None):
        if not paths:
            raise ValueError("At least one BAM file is required")
        if pysam is None:
            raise RuntimeError("pysam is required for BAM fragment input")
        self.input_paths = [str(Path(path)) for path in paths]
        require_bam_indexes(self.input_paths)
        self.handles = [pysam.AlignmentFile(path, "rb") for path in self.input_paths]
        try:
            merged = merge_bam_reference_headers_with_aliases(self.handles)
            supplied_sizes = _read_chrom_sizes(chrom_sizes)
            if supplied_sizes:
                supplied_names = list(supplied_sizes)
                seen_keys: dict[str, str] = {}
                for name in supplied_names:
                    key = canonical_contig_key(name)
                    previous = seen_keys.get(key)
                    if previous is not None and previous != name:
                        raise ValueError(
                            "Chromosome sizes contain ambiguous equivalent contigs "
                            f"{previous!r} and {name!r}."
                        )
                    seen_keys[key] = name

                merged_lengths = dict(zip(merged.references, merged.lengths))
                canonical_to_merged: dict[str, str] = {}
                for canonical in supplied_names:
                    try:
                        merged_name = resolve_contig_name(
                            canonical, merged.references, source_label="BAM headers"
                        )
                    except KeyError:
                        continue
                    if int(supplied_sizes[canonical]) != int(merged_lengths[merged_name]):
                        raise ValueError(
                            f"BAM contig {merged_name!r} has length "
                            f"{merged_lengths[merged_name]:,}, but chromosome sizes "
                            f"contig {canonical!r} has length {supplied_sizes[canonical]:,}."
                        )
                    canonical_to_merged[canonical] = merged_name

                self.references = supplied_names
                self.lengths = [int(supplied_sizes[name]) for name in supplied_names]
                self.source_contigs = []
                for source_mapping in merged.source_contigs:
                    self.source_contigs.append({
                        canonical: source_mapping[merged_name]
                        for canonical, merged_name in canonical_to_merged.items()
                        if merged_name in source_mapping
                    })
            else:
                self.references = merged.references
                self.lengths = merged.lengths
                self.source_contigs = merged.source_contigs
        except Exception:
            self.close()
            raise

    def fetch(
        self,
        contig: str,
        start: int,
        end: int,
        *,
        max_per_coordinate: int,
        subsample: float | None,
        dedup_scope: str,
    ) -> list[Fragment]:
        try:
            canonical = resolve_contig_name(contig, self.references, source_label="BAM headers")
        except KeyError as exc:
            raise ValueError(str(exc)) from exc
        source_names = [mapping.get(canonical) for mapping in self.source_contigs]
        return collect_fragments(
            bamfiles=self.handles,
            contig=canonical,
            start=start,
            end=end,
            max_per_coordinate=max_per_coordinate,
            subsample=subsample,
            dedup_scope=dedup_scope,
            source_contigs=source_names,
        )

    def close(self) -> None:
        for handle in getattr(self, "handles", []):
            try:
                handle.close()
            except Exception:
                pass


class IndexedIntervalFragmentSource(FragmentSource):
    """Random-access fragment source for tabix and bigBed inputs.

    Unlike :class:`IntervalFragmentSource`, this implementation does not first
    materialise every input record into SQLite.  Each worker opens the existing
    interval index and fetches only its assigned contig/range, which makes true
    per-contig multiprocessing possible without rescanning the whole file.
    """

    kind = "fragments"

    def __init__(
        self,
        paths: Sequence[str],
        *,
        chrom_sizes: str | None = None,
        fasta: pysam.FastaFile | None = None,
    ):
        if not paths:
            raise ValueError("At least one indexed fragment interval file is required")
        self.input_paths = [str(Path(path)) for path in paths]
        self._handles: list[tuple[str, object, list[str], dict[str, int]]] = []
        try:
            for raw_path in self.input_paths:
                path = Path(raw_path)
                if not path.exists():
                    raise FileNotFoundError(path)
                lower = path.name.lower()
                if lower.endswith((".bb", ".bigbed")):
                    if pyBigWig is None:
                        raise RuntimeError("pyBigWig is required for bigBed fragment input")
                    handle = pyBigWig.open(str(path))
                    sizes = {str(name): int(length) for name, length in handle.chroms().items()}
                    self._handles.append(("bigbed", handle, list(sizes), sizes))
                else:
                    if pysam is None:
                        raise RuntimeError("pysam is required for tabix-indexed fragment input")
                    if not (Path(str(path) + ".tbi").exists() or Path(str(path) + ".csi").exists()):
                        raise ValueError(f"Indexed interval input lacks .tbi/.csi index: {path}")
                    handle = pysam.TabixFile(str(path))
                    names = [str(name) for name in handle.contigs]
                    self._handles.append(("tabix", handle, names, {}))

            supplied_sizes = _read_chrom_sizes(chrom_sizes)
            fasta_sizes = (
                {str(name): int(length) for name, length in zip(fasta.references, fasta.lengths)}
                if fasta is not None
                else {}
            )
            if supplied_sizes:
                self.references = list(supplied_sizes)
                self.lengths = [int(supplied_sizes[name]) for name in self.references]
            elif fasta_sizes:
                available = {
                    canonical_contig_key(name)
                    for _kind, _handle, names, _sizes in self._handles
                    for name in names
                }
                self.references = [
                    name for name in fasta.references if canonical_contig_key(name) in available
                ]
                if not self.references:
                    raise ValueError("No indexed fragment contigs match the FASTA index")
                self.lengths = [int(fasta_sizes[name]) for name in self.references]
            else:
                # bigBed embeds chromosome lengths.  Tabix does not, so a
                # reference length source is required when any tabix input is used.
                if any(kind == "tabix" for kind, *_rest in self._handles):
                    raise ValueError(
                        "Tabix-indexed fragment input requires --chrom-sizes or --fasta "
                        "so chromosome lengths are known without scanning the file"
                    )
                names_by_key: dict[str, list[str]] = {}
                size_by_key: dict[str, int] = {}
                key_order: list[str] = []
                for _kind, _handle, names, sizes in self._handles:
                    for name in names:
                        key = canonical_contig_key(name)
                        if key not in names_by_key:
                            key_order.append(key)
                            names_by_key[key] = []
                        if name not in names_by_key[key]:
                            names_by_key[key].append(name)
                        length = int(sizes[name])
                        previous = size_by_key.get(key)
                        if previous is not None and previous != length:
                            raise ValueError(
                                f"Indexed fragment sources disagree on length for {name!r}: "
                                f"{previous:,} versus {length:,}"
                            )
                        size_by_key[key] = length
                self.references = [_choose_canonical_name(names_by_key[key]) for key in key_order]
                self.lengths = [size_by_key[key] for key in key_order]

            self._source_contigs: list[dict[str, str]] = []
            length_by_canonical = dict(zip(self.references, self.lengths))
            for kind, _handle, names, embedded_sizes in self._handles:
                mapping: dict[str, str] = {}
                for canonical in self.references:
                    try:
                        source_name = resolve_contig_name(
                            canonical, names, source_label="indexed fragment input"
                        )
                    except KeyError:
                        continue
                    if kind == "bigbed":
                        embedded = int(embedded_sizes[source_name])
                        expected = int(length_by_canonical[canonical])
                        if embedded != expected:
                            raise ValueError(
                                f"bigBed contig {source_name!r} has length {embedded:,}, "
                                f"but the analysis namespace gives {canonical!r} length {expected:,}"
                            )
                    mapping[canonical] = source_name
                self._source_contigs.append(mapping)
            if not any(self._source_contigs):
                raise ValueError("No indexed fragment contigs match the analysis namespace")
        except Exception:
            self.close()
            raise

    @staticmethod
    def _parse_interval(raw: str, path: str) -> Fragment:
        fields = raw.rstrip("\r\n").split("\t")
        if len(fields) < 3:
            fields = raw.split()
        if len(fields) < 3:
            raise ValueError(f"{path}: indexed record has fewer than three columns")
        try:
            start = int(fields[1]); end = int(fields[2])
        except ValueError as exc:
            raise ValueError(f"{path}: indexed start/end must be integers") from exc
        if start < 0 or end <= start:
            raise ValueError(f"{path}: indexed record requires 0 <= start < end")
        return start, end

    def fetch(
        self,
        contig: str,
        start: int,
        end: int,
        *,
        max_per_coordinate: int,
        subsample: float | None,
        dedup_scope: str,
    ) -> list[Fragment]:
        if max_per_coordinate < 0:
            raise ValueError("max_per_coordinate must be 0 or greater")
        if dedup_scope not in {"all_bams", "per_bam"}:
            raise ValueError("dedup_scope must be 'all_bams' or 'per_bam'")
        if subsample is not None and not 0.0 <= subsample <= 1.0:
            raise ValueError("subsample must be between 0 and 1")
        try:
            canonical = resolve_contig_name(
                contig, self.references, source_label="indexed fragment inputs"
            )
        except KeyError as exc:
            raise ValueError(str(exc)) from exc

        records: list[tuple[int, int, int]] = []
        for source_id, ((kind, handle, _names, _sizes), mapping, path) in enumerate(
            zip(self._handles, self._source_contigs, self.input_paths)
        ):
            source_name = mapping.get(canonical)
            if source_name is None:
                continue
            if kind == "tabix":
                try:
                    iterator = handle.fetch(source_name, max(0, start), max(start, end))
                except ValueError as exc:
                    raise ValueError(f"Unable to fetch {source_name!r} from {path}: {exc}") from exc
                for raw in iterator:
                    fragment_start, fragment_end = self._parse_interval(raw, path)
                    if fragment_start < end and fragment_end > start:
                        records.append((source_id, fragment_start, fragment_end))
            else:
                entries = handle.entries(source_name, max(0, start), max(start, end)) or []
                for fragment_start, fragment_end, _rest in entries:
                    if int(fragment_start) < end and int(fragment_end) > start:
                        records.append((source_id, int(fragment_start), int(fragment_end)))

        records.sort(key=lambda item: (item[1], item[2], item[0]))
        counts: dict[tuple[int, ...], int] = defaultdict(int)
        output: list[Fragment] = []
        for source_id, fragment_start, fragment_end in records:
            key = (
                (fragment_start, fragment_end)
                if dedup_scope == "all_bams"
                else (source_id, fragment_start, fragment_end)
            )
            if max_per_coordinate > 0 and counts[key] >= max_per_coordinate:
                continue
            counts[key] += 1
            if subsample is not None and random.random() > subsample:
                continue
            output.append((fragment_start, fragment_end))
        return output

    def iter_all(
        self,
        *,
        contigs: set[str] | None = None,
        min_length: int = 1,
        max_length: int | None = None,
        max_per_coordinate: int = 0,
        dedup_scope: str = "all_bams",
        subsample: float | None = None,
    ) -> Iterator[tuple[str, int, int]]:
        selected = list(self.references) if contigs is None else []
        if contigs is not None:
            for requested in contigs:
                try:
                    canonical = resolve_contig_name(
                        requested, self.references, source_label="indexed fragment inputs"
                    )
                except KeyError as exc:
                    raise ValueError(str(exc)) from exc
                if canonical not in selected:
                    selected.append(canonical)
        length_by_name = dict(zip(self.references, self.lengths))
        for chrom in selected:
            for start, end in self.fetch(
                chrom, 0, length_by_name[chrom],
                max_per_coordinate=max_per_coordinate,
                subsample=subsample,
                dedup_scope=dedup_scope,
            ):
                size = end - start
                if size < min_length or (max_length is not None and size > max_length):
                    continue
                yield chrom, start, end

    def close(self) -> None:
        for _kind, handle, _names, _sizes in getattr(self, "_handles", []):
            try:
                handle.close()
            except Exception:
                pass
        self._handles = []


class IntervalFragmentSource(FragmentSource):
    """Disk-indexed BED/BED.gz/bigBed fragment source.

    The SQLite table retains both the original source spelling and the canonical
    analysis spelling.  Queries and duplicate counting use canonical names, so
    records labelled ``20`` and ``chr20`` can participate in one analysis when
    the supplied chromosome-size namespace identifies them as equivalent.
    """

    kind = "fragments"

    def __init__(
        self,
        paths: Sequence[str],
        *,
        chrom_sizes: str | None = None,
        fasta: pysam.FastaFile | None = None,
    ):
        if not paths:
            raise ValueError("At least one fragment BED, BED.gz or bigBed file is required")
        self.input_paths = [str(Path(path)) for path in paths]
        for path in self.input_paths:
            if not Path(path).exists():
                raise FileNotFoundError(path)

        supplied_sizes = _read_chrom_sizes(chrom_sizes)
        supplied_names = list(supplied_sizes)
        fasta_sizes = dict(zip(fasta.references, fasta.lengths)) if fasta is not None else {}

        temp = tempfile.NamedTemporaryFile(
            prefix="nucleosuite_fragments_", suffix=".sqlite", delete=False
        )
        temp.close()
        self._db_path = temp.name
        self._connection = sqlite3.connect(self._db_path)
        self._connection.execute("PRAGMA journal_mode=OFF")
        self._connection.execute("PRAGMA synchronous=OFF")
        self._connection.execute("PRAGMA temp_store=MEMORY")
        self._connection.execute(
            "CREATE TABLE fragments ("
            "source_id INTEGER NOT NULL, source_chrom TEXT NOT NULL, "
            "canonical_chrom TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL)"
        )
        self._connection.execute(
            "CREATE TABLE contig_aliases ("
            "source_id INTEGER NOT NULL, source_chrom TEXT NOT NULL, "
            "canonical_chrom TEXT NOT NULL, length INTEGER, "
            "PRIMARY KEY (source_id, source_chrom))"
        )

        key_order: list[str] = []
        names_by_key: dict[str, list[str]] = {}
        max_end_by_name: dict[str, int] = {}
        source_name_by_key: list[dict[str, str]] = []
        canonical_for_name: dict[str, str] = {}
        if supplied_sizes:
            for name in supplied_names:
                canonical_for_name[name] = name

        batch: list[tuple[int, str, str, int, int]] = []
        valid_count = 0
        try:
            for source_id, path in enumerate(self.input_paths):
                current: dict[str, str] = {}
                with open_text(path) as handle:
                    for line_number, raw in enumerate(handle, 1):
                        text = raw.strip()
                        if not text or text.startswith(("#", "track", "browser")):
                            continue
                        fields = text.split("\t") if "\t" in text else text.split()
                        if len(fields) < 3:
                            raise ValueError(f"{path}:{line_number}: expected at least three columns")
                        source_chrom = fields[0]
                        try:
                            start = int(fields[1])
                            end = int(fields[2])
                        except ValueError as exc:
                            raise ValueError(
                                f"{path}:{line_number}: start and end must be integers"
                            ) from exc
                        if start < 0 or end <= start:
                            raise ValueError(f"{path}:{line_number}: require 0 <= start < end")

                        key = canonical_contig_key(source_chrom)
                        previous = current.get(key)
                        if previous is not None and previous != source_chrom:
                            raise ValueError(
                                f"{path} contains ambiguous equivalent contigs "
                                f"{previous!r} and {source_chrom!r}."
                            )
                        current[key] = source_chrom
                        if key not in names_by_key:
                            key_order.append(key)
                            names_by_key[key] = [source_chrom]
                        elif source_chrom not in names_by_key[key]:
                            names_by_key[key].append(source_chrom)
                        max_end_by_name[source_chrom] = max(
                            max_end_by_name.get(source_chrom, 0), end
                        )

                        if supplied_sizes:
                            try:
                                canonical = resolve_contig_name(
                                    source_chrom,
                                    supplied_names,
                                    source_label="chromosome sizes",
                                )
                            except KeyError as exc:
                                raise ValueError(
                                    f"{path}:{line_number}: {str(exc)}"
                                ) from exc
                        else:
                            canonical = source_chrom  # normalized after discovery
                        batch.append((source_id, source_chrom, canonical, start, end))
                        valid_count += 1
                        if len(batch) >= 100_000:
                            self._connection.executemany(
                                "INSERT INTO fragments VALUES (?,?,?,?,?)", batch
                            )
                            batch.clear()
                source_name_by_key.append(current)
            if batch:
                self._connection.executemany("INSERT INTO fragments VALUES (?,?,?,?,?)", batch)
            if valid_count == 0:
                raise ValueError("No valid fragment intervals were found")

            if supplied_sizes:
                self.references = list(supplied_sizes)
                self.lengths = [int(supplied_sizes[name]) for name in self.references]
                for source_id, current in enumerate(source_name_by_key):
                    for source_name in current.values():
                        canonical = resolve_contig_name(
                            source_name, self.references, source_label="chromosome sizes"
                        )
                        canonical_for_name[source_name] = canonical
                        maximum = max_end_by_name[source_name]
                        if maximum > supplied_sizes[canonical]:
                            raise ValueError(
                                f"Fragment interval on {source_name} ends at {maximum:,}, "
                                f"beyond supplied chromosome length {supplied_sizes[canonical]:,}"
                            )
                        self._connection.execute(
                            "INSERT OR REPLACE INTO contig_aliases VALUES (?,?,?,?)",
                            (source_id, source_name, canonical, int(supplied_sizes[canonical])),
                        )
            else:
                canonical_by_key = {
                    key: _choose_canonical_name(names_by_key[key]) for key in key_order
                }
                self.references = [canonical_by_key[key] for key in key_order]
                lengths: list[int] = []
                for key in key_order:
                    canonical = canonical_by_key[key]
                    maximum = max(max_end_by_name[name] for name in names_by_key[key])
                    length = None
                    if fasta_sizes:
                        try:
                            fasta_name = resolve_contig_name(
                                canonical, list(fasta_sizes), source_label="FASTA index"
                            )
                            length = int(fasta_sizes[fasta_name])
                        except KeyError:
                            length = None
                    if length is None:
                        length = maximum
                    if length < maximum:
                        raise ValueError(
                            f"Fragment interval on {canonical} ends at {maximum:,}, "
                            f"beyond chromosome length {length:,}"
                        )
                    lengths.append(length)
                    for source_name in names_by_key[key]:
                        canonical_for_name[source_name] = canonical
                self.lengths = lengths

                for source_id, current in enumerate(source_name_by_key):
                    for key, source_name in current.items():
                        canonical = canonical_by_key[key]
                        length = self.lengths[self.references.index(canonical)]
                        self._connection.execute(
                            "UPDATE fragments SET canonical_chrom=? "
                            "WHERE source_id=? AND source_chrom=?",
                            (canonical, source_id, source_name),
                        )
                        self._connection.execute(
                            "INSERT OR REPLACE INTO contig_aliases VALUES (?,?,?,?)",
                            (source_id, source_name, canonical, int(length)),
                        )

            self._connection.execute(
                "CREATE INDEX fragment_lookup ON fragments(canonical_chrom, start, end)"
            )
            self._connection.commit()
        except Exception:
            self.close()
            raise

    def fetch(
        self,
        contig: str,
        start: int,
        end: int,
        *,
        max_per_coordinate: int,
        subsample: float | None,
        dedup_scope: str,
    ) -> list[Fragment]:
        if max_per_coordinate < 0:
            raise ValueError("max_per_coordinate must be 0 or greater")
        if dedup_scope not in {"all_bams", "per_bam"}:
            raise ValueError("dedup_scope must be 'all_bams' or 'per_bam'")
        if subsample is not None and not 0.0 <= subsample <= 1.0:
            raise ValueError("subsample must be between 0 and 1")
        try:
            canonical = resolve_contig_name(
                contig, self.references, source_label="fragment inputs"
            )
        except KeyError as exc:
            raise ValueError(str(exc)) from exc

        cursor = self._connection.execute(
            "SELECT source_id, start, end FROM fragments "
            "WHERE canonical_chrom=? AND start<? AND end>? "
            "ORDER BY start, end, source_id",
            (canonical, end, start),
        )
        counts: dict[tuple[int, ...], int] = defaultdict(int)
        output: list[Fragment] = []
        for source_id, fragment_start, fragment_end in cursor:
            key = (
                (int(fragment_start), int(fragment_end))
                if dedup_scope == "all_bams"
                else (int(source_id), int(fragment_start), int(fragment_end))
            )
            if max_per_coordinate > 0 and counts[key] >= max_per_coordinate:
                continue
            counts[key] += 1
            if subsample is not None and random.random() > subsample:
                continue
            output.append((int(fragment_start), int(fragment_end)))
        return output

    def iter_all(
        self,
        *,
        contigs: set[str] | None = None,
        min_length: int = 1,
        max_length: int | None = None,
        max_per_coordinate: int = 0,
        dedup_scope: str = "all_bams",
        subsample: float | None = None,
    ) -> Iterator[tuple[str, int, int]]:
        if contigs is None:
            selected = list(self.references)
        else:
            selected = []
            for requested in contigs:
                try:
                    canonical = resolve_contig_name(
                        requested, self.references, source_label="fragment inputs"
                    )
                except KeyError as exc:
                    raise ValueError(str(exc)) from exc
                if canonical not in selected:
                    selected.append(canonical)
        length_by_name = dict(zip(self.references, self.lengths))
        for chrom in selected:
            for start, end in self.fetch(
                chrom,
                0,
                length_by_name[chrom],
                max_per_coordinate=max_per_coordinate,
                subsample=subsample,
                dedup_scope=dedup_scope,
            ):
                size = end - start
                if size < min_length or (max_length is not None and size > max_length):
                    continue
                yield chrom, start, end

    def close(self) -> None:
        connection = getattr(self, "_connection", None)
        if connection is not None:
            try:
                connection.close()
            except Exception:
                pass
            self._connection = None
        path = getattr(self, "_db_path", None)
        if path:
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass
            self._db_path = None


def open_fragment_source(
    *,
    bam_paths: Sequence[str] | None = None,
    fragment_paths: Sequence[str] | None = None,
    chrom_sizes: str | None = None,
    fasta: pysam.FastaFile | None = None,
    blacklist_path: str | None = None,
) -> FragmentSource:
    if bool(bam_paths) == bool(fragment_paths):
        raise ValueError("Provide exactly one of --bamfiles/--bam or --fragments")
    if bam_paths:
        source: FragmentSource = BamFragmentSource(bam_paths, chrom_sizes=chrom_sizes)
    else:
        paths = [Path(value) for value in (fragment_paths or [])]
        indexed = bool(paths) and all(
            path.name.lower().endswith((".bb", ".bigbed"))
            or (
                path.name.lower().endswith((".bed.gz", ".bed.bgz", ".tsv.gz", ".tsv.bgz"))
                and (Path(str(path) + ".tbi").exists() or Path(str(path) + ".csi").exists())
            )
            for path in paths
        )
        if indexed:
            source = IndexedIntervalFragmentSource(
                [str(path) for path in paths], chrom_sizes=chrom_sizes, fasta=fasta
            )
        else:
            source = IntervalFragmentSource(
                [str(path) for path in paths], chrom_sizes=chrom_sizes, fasta=fasta
            )
    try:
        blacklist = load_blacklist(blacklist_path, source.references, source.lengths)
        if blacklist is not None:
            source = BlacklistFragmentSource(source, blacklist)
        return source
    except Exception:
        source.close()
        raise
