"""Create mutually exclusive gene categories from gene and chromatin-state BED files."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.core.chrom_sizes import chromosome_size_dict
from nucleosuite.io import open_text
from nucleosuite.io.intervals import INTERVAL_FORMATS, finalise_interval_files
from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.core.blacklist import BlacklistIndex, load_blacklist_unbounded
from nucleosuite.parallel import add_parallel_arguments
from nucleosuite.partitioned import run_partitioned_command
from nucleosuite.progress import ProgressReporter


_TOKEN_RE = re.compile(r"\s*([&|()]|[^&|()\s]+)")
_SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9._-]+")


@dataclass(frozen=True)
class GeneRecord:
    chrom: str
    start: int
    end: int
    gene_id: str
    fields: tuple[str, ...]


@dataclass(frozen=True)
class StateRecord:
    chrom: str
    start: int
    end: int
    label: str


@dataclass(frozen=True)
class GeneSetRule:
    name: str
    expression: str
    rpn: tuple[str, ...]
    state_names: frozenset[str]
    exclude_if_candidate: frozenset[str] = frozenset()


def safe_name(value: str) -> str:
    cleaned = _SAFE_NAME_RE.sub("_", value.strip()).strip("_")
    if not cleaned:
        raise ValueError(f"Unable to make a safe file name from {value!r}")
    return cleaned


def _iter_bed_fields(path: str | Path) -> Iterable[tuple[int, list[str]]]:
    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t") if "\t" in text else text.split()
            yield line_no, fields


def read_chrom_sizes(path: str | Path | None) -> dict[str, int] | None:
    if path is None:
        return None
    return chromosome_size_dict(path)


def read_genes(
    path: str | Path,
    gene_id_column: int = 4,
    chrom_sizes: dict[str, int] | None = None,
) -> list[GeneRecord]:
    if gene_id_column < 1:
        raise ValueError("gene ID column must be one-based and positive")
    idx = gene_id_column - 1
    genes: list[GeneRecord] = []
    seen_ids: dict[str, tuple[str, int, int]] = {}
    for line_no, fields in _iter_bed_fields(path):
        if len(fields) < 3:
            raise ValueError(f"{path}:{line_no}: expected BED3 or greater")
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError as exc:
            raise ValueError(f"{path}:{line_no}: invalid BED coordinates") from exc
        chrom = fields[0]
        if start < 0 or end <= start:
            raise ValueError(f"{path}:{line_no}: require 0 <= start < end")
        if chrom_sizes is not None:
            try:
                chrom = resolve_contig_name(
                    chrom, list(chrom_sizes), source_label="chromosome sizes"
                )
            except KeyError:
                continue
            end = min(end, chrom_sizes[chrom])
            if end <= start:
                continue
        gene_id = fields[idx] if idx < len(fields) and fields[idx] else f"{chrom}:{start}-{end}"
        location = (chrom, start, end)
        previous = seen_ids.get(gene_id)
        if previous is not None:
            raise ValueError(
                f"Gene ID {gene_id!r} occurs more than once ({previous} and {location}). "
                "Use a gene BED with one unique record per gene."
            )
        seen_ids[gene_id] = location
        out_fields = list(fields)
        out_fields[0] = chrom
        out_fields[1] = str(start)
        out_fields[2] = str(end)
        if len(out_fields) < 4:
            out_fields.append(gene_id)
        genes.append(GeneRecord(chrom, start, end, gene_id, tuple(out_fields)))
    genes.sort(key=lambda g: (g.chrom, g.start, g.end, g.gene_id))
    if not genes:
        raise ValueError(f"No compatible genes found in {path}")
    return genes


def read_states(
    path: str | Path,
    state_label_column: int = 4,
    chrom_sizes: dict[str, int] | None = None,
) -> list[StateRecord]:
    if state_label_column < 1:
        raise ValueError("state label column must be one-based and positive")
    idx = state_label_column - 1
    states: list[StateRecord] = []
    for line_no, fields in _iter_bed_fields(path):
        if len(fields) <= idx or len(fields) < 3:
            raise ValueError(
                f"{path}:{line_no}: state label column {state_label_column} is unavailable"
            )
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError as exc:
            raise ValueError(f"{path}:{line_no}: invalid BED coordinates") from exc
        chrom = fields[0]
        if start < 0 or end <= start:
            raise ValueError(f"{path}:{line_no}: require 0 <= start < end")
        if chrom_sizes is not None:
            try:
                chrom = resolve_contig_name(
                    chrom, list(chrom_sizes), source_label="chromosome sizes"
                )
            except KeyError:
                continue
            end = min(end, chrom_sizes[chrom])
            if end <= start:
                continue
        states.append(StateRecord(chrom, start, end, fields[idx]))
    states.sort(key=lambda s: (s.chrom, s.start, s.end, s.label))
    if not states:
        raise ValueError(f"No compatible chromatin-state records found in {path}")
    return states


def tokenize(expression: str) -> list[str]:
    if "!" in expression:
        raise ValueError(
            "Gene-set rules support inclusion only and do not use '!'. Put candidate-set exclusions in "
            "the exclude_if_candidate configuration column."
        )
    tokens: list[str] = []
    pos = 0
    while pos < len(expression):
        match = _TOKEN_RE.match(expression, pos)
        if match is None:
            raise ValueError(f"Unable to parse rule near: {expression[pos:]}")
        tokens.append(match.group(1))
        pos = match.end()
    if not tokens:
        raise ValueError("Gene-set rule is empty")
    return tokens


def expression_to_rpn(expression: str) -> tuple[tuple[str, ...], frozenset[str]]:
    tokens = tokenize(expression)
    precedence = {"|": 1, "&": 2}
    output: list[str] = []
    operators: list[str] = []
    state_names: set[str] = set()
    expect_operand = True
    for token in tokens:
        if token in ("&", "|"):
            if expect_operand:
                raise ValueError(f"Unexpected operator {token!r} in rule {expression!r}")
            while operators and operators[-1] != "(" and precedence[operators[-1]] >= precedence[token]:
                output.append(operators.pop())
            operators.append(token)
            expect_operand = True
        elif token == "(":
            if not expect_operand:
                raise ValueError(f"Missing operator before '(' in rule {expression!r}")
            operators.append(token)
        elif token == ")":
            if expect_operand:
                raise ValueError(f"Unexpected ')' in rule {expression!r}")
            while operators and operators[-1] != "(":
                output.append(operators.pop())
            if not operators:
                raise ValueError(f"Unbalanced ')' in rule {expression!r}")
            operators.pop()
            expect_operand = False
        else:
            if not expect_operand:
                raise ValueError(f"Missing operator before {token!r} in rule {expression!r}")
            output.append(token)
            state_names.add(token)
            expect_operand = False
    if expect_operand:
        raise ValueError(f"Rule ends with an operator: {expression!r}")
    while operators:
        operator = operators.pop()
        if operator == "(":
            raise ValueError(f"Unbalanced '(' in rule {expression!r}")
        output.append(operator)
    return tuple(output), frozenset(state_names)


def evaluate_rpn(rpn: Sequence[str], intersecting_states: set[str]) -> bool:
    stack: list[bool] = []
    for token in rpn:
        if token in ("&", "|"):
            if len(stack) < 2:
                raise ValueError("Invalid rule expression")
            right = stack.pop()
            left = stack.pop()
            stack.append(left and right if token == "&" else left or right)
        else:
            stack.append(token in intersecting_states)
    if len(stack) != 1:
        raise ValueError("Invalid rule expression")
    return stack[0]


def _split_set_names(value: str | None) -> frozenset[str]:
    if not value:
        return frozenset()
    return frozenset(
        item.strip()
        for item in re.split(r"[,;]", value)
        if item.strip()
    )


def load_rules(
    config_path: str | Path | None,
    inline_rules: Sequence[str] | None,
) -> list[GeneSetRule]:
    rows: list[tuple[str, str, frozenset[str]]] = []
    if config_path is not None and inline_rules:
        raise ValueError("Use either --config or --gene-set, not both")
    if inline_rules:
        for item in inline_rules:
            if "=" not in item:
                raise ValueError("--gene-set values must use NAME=RULE")
            name, expression = item.split("=", 1)
            rows.append((name.strip(), expression.strip(), frozenset()))
    elif config_path is not None:
        with open(config_path, encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or not {"set_name", "include_rule"} <= set(reader.fieldnames):
                raise ValueError("Gene-set config requires columns: set_name and include_rule")
            exclusion_column = next(
                (
                    name
                    for name in (
                        "exclude_if_candidate",
                        "exclude_candidate_sets",
                        "exclude_if_candidates",
                    )
                    if name in reader.fieldnames
                ),
                None,
            )
            for row in reader:
                name = (row.get("set_name") or "").strip()
                expression = (row.get("include_rule") or "").strip()
                exclusions = _split_set_names(row.get(exclusion_column) if exclusion_column else None)
                if name or expression or exclusions:
                    rows.append((name, expression, exclusions))
    else:
        raise ValueError("A gene-set config or one or more --gene-set rules is required")

    rules: list[GeneSetRule] = []
    seen: set[str] = set()
    for name, expression, exclusions in rows:
        if not name or not expression:
            raise ValueError("Each gene-set definition requires a name and include rule")
        if name in seen:
            raise ValueError(f"Duplicate gene-set name: {name}")
        seen.add(name)
        rpn, state_names = expression_to_rpn(expression)
        rules.append(GeneSetRule(name, expression, rpn, state_names, exclusions))
    if len(rules) < 2:
        raise ValueError("Define at least two gene sets so overlap handling is meaningful")

    names = {rule.name for rule in rules}
    unknown = sorted(
        excluded
        for rule in rules
        for excluded in rule.exclude_if_candidate
        if excluded not in names
    )
    if unknown:
        raise ValueError(
            "Unknown candidate set(s) in exclude_if_candidate: " + ", ".join(sorted(set(unknown)))
        )
    for rule in rules:
        if rule.name in rule.exclude_if_candidate:
            raise ValueError(f"Set {rule.name!r} cannot exclude itself")
    return rules


def intersect_states_by_gene(
    genes: Sequence[GeneRecord], states: Sequence[StateRecord]
) -> dict[str, set[str]]:
    genes_by_chrom: dict[str, list[GeneRecord]] = defaultdict(list)
    states_by_chrom: dict[str, list[StateRecord]] = defaultdict(list)
    for gene in genes:
        genes_by_chrom[gene.chrom].append(gene)
    for state in states:
        states_by_chrom[state.chrom].append(state)

    result: dict[str, set[str]] = {gene.gene_id: set() for gene in genes}
    for chrom, chrom_genes in genes_by_chrom.items():
        chrom_states = states_by_chrom.get(chrom, [])
        state_index = 0
        active: list[StateRecord] = []
        for gene in chrom_genes:
            active = [state for state in active if state.end > gene.start]
            while state_index < len(chrom_states) and chrom_states[state_index].start < gene.end:
                state = chrom_states[state_index]
                if state.end > gene.start:
                    active.append(state)
                state_index += 1
            result[gene.gene_id].update(
                state.label
                for state in active
                if state.start < gene.end and state.end > gene.start
            )
    return result


def write_bed(path: Path, records: Sequence[GeneRecord]) -> None:
    """Write a standards-compliant BED6 gene interval file.

    The source gene BED may contain arbitrary annotations after column 3. Those
    annotations are retained in the gene-assignment TSV, while interval outputs
    use BED6 so they can always be converted to bigBed without misinterpreting a
    gene name as the numeric BED score field.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for record in records:
            handle.write(
                "\t".join(
                    [
                        record.chrom,
                        str(record.start),
                        str(record.end),
                        record.gene_id,
                        "0",
                        _gene_strand(record),
                    ]
                )
                + "\n"
            )


def write_tss_bed(path: Path, records: Sequence[GeneRecord]) -> None:
    """Write one-base BED6 transcription-start-site intervals.

    Plus-strand genes use ``start`` to ``start + 1``. Minus-strand genes use
    ``end - 1`` to ``end``. Records without a valid strand are omitted because
    they cannot be oriented for strand-aware aggregation.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for record in records:
            strand = _gene_strand(record)
            if strand == "+":
                start, end = record.start, record.start + 1
            elif strand == "-":
                start, end = record.end - 1, record.end
            else:
                continue
            handle.write(
                "\t".join(
                    [record.chrom, str(start), str(end), record.gene_id, "0", strand]
                )
                + "\n"
            )


def _gene_name(record: GeneRecord) -> str:
    return record.fields[4] if len(record.fields) >= 5 and record.fields[4] else record.gene_id


def _gene_strand(record: GeneRecord) -> str:
    return record.fields[5] if len(record.fields) >= 6 and record.fields[5] in {"+", "-"} else "."


def gene_anchor_interval(record: GeneRecord) -> tuple[int, int]:
    """Return the one-base TSS anchor, or the full interval without a strand."""
    strand = _gene_strand(record)
    if strand == "+":
        return record.start, record.start + 1
    if strand == "-":
        return record.end - 1, record.end
    return record.start, record.end


def filter_blacklisted_gene_anchors(
    genes: Sequence[GeneRecord], blacklist: BlacklistIndex | None
) -> tuple[list[GeneRecord], int]:
    """Exclude genes whose TSS anchor overlaps a blacklist interval."""
    if blacklist is None:
        return list(genes), 0
    retained: list[GeneRecord] = []
    excluded = 0
    for gene in genes:
        start, end = gene_anchor_interval(gene)
        if blacklist.overlaps(gene.chrom, start, end):
            excluded += 1
        else:
            retained.append(gene)
    return retained, excluded


def write_state_labeled_bed(
    path: Path,
    genes: Sequence[GeneRecord],
    final_ids: dict[str, set[str]],
    set_order: Sequence[str],
) -> None:
    """Write BED6 with the final category in column 4 for state-aware DAC."""
    category_by_gene: dict[str, str] = {}
    for set_name in set_order:
        for gene_id in final_ids.get(set_name, set()):
            if gene_id in category_by_gene:
                raise ValueError(
                    f"Gene {gene_id!r} remains assigned to more than one final set: "
                    f"{category_by_gene[gene_id]!r} and {set_name!r}"
                )
            category_by_gene[gene_id] = set_name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for gene in genes:
            category = category_by_gene.get(gene.gene_id)
            if category is None:
                continue
            handle.write(
                "\t".join(
                    [
                        gene.chrom,
                        str(gene.start),
                        str(gene.end),
                        category,
                        "0",
                        _gene_strand(gene),
                    ]
                )
                + "\n"
            )


def write_state_labeled_tss_bed(
    path: Path,
    genes: Sequence[GeneRecord],
    final_ids: dict[str, set[str]],
    set_order: Sequence[str],
) -> None:
    """Write BED6 TSS intervals with the final gene-set label in column 4."""
    category_by_gene: dict[str, str] = {}
    for set_name in set_order:
        for gene_id in final_ids.get(set_name, set()):
            if gene_id in category_by_gene:
                raise ValueError(
                    f"Gene {gene_id!r} remains assigned to more than one final set: "
                    f"{category_by_gene[gene_id]!r} and {set_name!r}"
                )
            category_by_gene[gene_id] = set_name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for gene in genes:
            category = category_by_gene.get(gene.gene_id)
            if category is None:
                continue
            strand = _gene_strand(gene)
            if strand == "+":
                start, end = gene.start, gene.start + 1
            elif strand == "-":
                start, end = gene.end - 1, gene.end
            else:
                continue
            handle.write(
                "\t".join(
                    [gene.chrom, str(start), str(end), category, "0", strand]
                )
                + "\n"
            )


def create_venn(
    set_names: Sequence[str],
    candidate_ids: dict[str, set[str]],
    output_path: Path,
    title: str,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from matplotlib_venn import venn2, venn3
    except ImportError as exc:
        raise RuntimeError(
            "Venn-diagram output requires matplotlib-venn. Update the NucleoSuite "
            "environment with 'mamba env update -n nucleosuite -f environment.yml'."
        ) from exc

    selected = [candidate_ids[name] for name in set_names]
    figure = plt.figure(figsize=(8, 8))
    if len(selected) == 2:
        venn2(selected, set_labels=list(set_names))
    elif len(selected) == 3:
        venn3(selected, set_labels=list(set_names))
    else:
        raise ValueError("Venn diagrams require exactly two or three selected gene sets")
    plt.title(title)
    plt.tight_layout()
    from nucleosuite.plotting import save_figure
    saved = save_figure(figure, output_path, default_dpi=300, bbox_inches="tight")
    plt.close(figure)
    return saved


def build_gene_sets(
    genes: Sequence[GeneRecord],
    states: Sequence[StateRecord],
    rules: Sequence[GeneSetRule],
    output_dir: str | Path,
    output_prefix: str = "gene_sets",
    venn_sets: Sequence[str] | None = None,
    interval_format: str = "bed",
    chrom_sizes=None,
    leftover_set_name: str | None = None,
    prefix_member_files: bool = False,
) -> dict[str, Path]:
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = safe_name(output_prefix)
    candidate_dir = outdir / "candidate_sets"
    final_dir = outdir / "final_sets"
    final_tss_dir = outdir / "final_tss"
    overlap_dir = outdir / "overlaps"
    for directory in (candidate_dir, final_dir, final_tss_dir, overlap_dir):
        directory.mkdir(parents=True, exist_ok=True)

    state_by_gene = intersect_states_by_gene(genes, states)
    available_states = {state.label for state in states}
    missing_states = sorted(
        {state_name for rule in rules for state_name in rule.state_names} - available_states
    )
    if missing_states:
        print(
            "WARNING: gene-set rules reference state labels absent from the compatible "
            "state records; those labels evaluate as false: " + ", ".join(missing_states),
            file=sys.stderr,
        )

    candidate_ids: dict[str, set[str]] = {}
    for rule in rules:
        candidate_ids[rule.name] = {
            gene.gene_id
            for gene in genes
            if evaluate_rpn(rule.rpn, state_by_gene[gene.gene_id])
        }

    membership: dict[str, list[str]] = defaultdict(list)
    for set_name, identifiers in candidate_ids.items():
        for gene_id in identifiers:
            membership[gene_id].append(set_name)
    overlapping_ids = {gene_id for gene_id, names in membership.items() if len(names) > 1}

    directed_exclusions = any(rule.exclude_if_candidate for rule in rules)
    excluded_ids_by_set: dict[str, set[str]] = {}
    if directed_exclusions:
        final_ids: dict[str, set[str]] = {}
        for rule in rules:
            excluded = set().union(
                *(candidate_ids[name] for name in rule.exclude_if_candidate)
            ) if rule.exclude_if_candidate else set()
            excluded_ids_by_set[rule.name] = candidate_ids[rule.name] & excluded
            final_ids[rule.name] = candidate_ids[rule.name] - excluded
    else:
        final_ids = {
            set_name: identifiers - overlapping_ids
            for set_name, identifiers in candidate_ids.items()
        }
        excluded_ids_by_set = {
            set_name: identifiers & overlapping_ids
            for set_name, identifiers in candidate_ids.items()
        }

    final_membership: dict[str, list[str]] = defaultdict(list)
    for set_name, identifiers in final_ids.items():
        for gene_id in identifiers:
            final_membership[gene_id].append(set_name)
    nonexclusive = {gene_id: names for gene_id, names in final_membership.items() if len(names) > 1}
    if nonexclusive:
        examples = "; ".join(
            f"{gene_id}={','.join(names)}" for gene_id, names in list(nonexclusive.items())[:5]
        )
        raise ValueError(
            "Final gene sets are not mutually exclusive. Add exclude_if_candidate rules "
            f"for the remaining overlaps. Examples: {examples}"
        )

    ordered_set_names = [rule.name for rule in rules]
    if leftover_set_name:
        leftover_set_name = safe_name(leftover_set_name)
        if leftover_set_name in final_ids:
            raise ValueError(f"Leftover set name duplicates a configured set: {leftover_set_name}")
        candidate_union = set().union(*candidate_ids.values()) if candidate_ids else set()
        final_ids[leftover_set_name] = {
            gene.gene_id for gene in genes
        } - candidate_union
        excluded_ids_by_set[leftover_set_name] = set()
        ordered_set_names.append(leftover_set_name)

    candidate_paths: dict[str, Path] = {}
    final_paths: dict[str, Path] = {}
    final_tss_paths: dict[str, Path] = {}
    member_prefix = f"{prefix}_" if prefix_member_files else ""
    for rule in rules:
        filename = f"{member_prefix}{safe_name(rule.name)}.bed"
        candidate_path = candidate_dir / filename
        final_path = final_dir / filename
        final_tss_path = final_tss_dir / filename
        selected_candidates = [gene for gene in genes if gene.gene_id in candidate_ids[rule.name]]
        selected_final = [gene for gene in genes if gene.gene_id in final_ids[rule.name]]
        write_bed(candidate_path, selected_candidates)
        write_bed(final_path, selected_final)
        write_tss_bed(final_tss_path, selected_final)
        candidate_paths[rule.name] = candidate_path
        final_paths[rule.name] = final_path
        final_tss_paths[rule.name] = final_tss_path
    if leftover_set_name:
        leftover_path = final_dir / f"{member_prefix}{safe_name(leftover_set_name)}.bed"
        leftover_tss_path = final_tss_dir / f"{member_prefix}{safe_name(leftover_set_name)}.bed"
        selected_leftover = [gene for gene in genes if gene.gene_id in final_ids[leftover_set_name]]
        write_bed(leftover_path, selected_leftover)
        write_tss_bed(leftover_tss_path, selected_leftover)
        final_paths[leftover_set_name] = leftover_path
        final_tss_paths[leftover_set_name] = leftover_tss_path

    overlap_path = overlap_dir / f"{member_prefix}shared_by_multiple_sets.bed"
    write_bed(overlap_path, [gene for gene in genes if gene.gene_id in overlapping_ids])

    final_state_path = outdir / f"{prefix}_final_states.bed"
    write_state_labeled_bed(final_state_path, genes, final_ids, ordered_set_names)
    final_tss_state_path = outdir / f"{prefix}_final_tss.bed"
    write_state_labeled_tss_bed(final_tss_state_path, genes, final_ids, ordered_set_names)

    interval_beds = [
        *candidate_paths.values(),
        *final_paths.values(),
        *final_tss_paths.values(),
        overlap_path,
        final_state_path,
        final_tss_state_path,
    ]
    if interval_format != "bed" and chrom_sizes is None:
        raise ValueError("--chrom-sizes is required for gene-set bigBed output")
    finalise_interval_files(interval_beds, interval_format, chrom_sizes or {})

    def selected_interval(path: Path) -> Path:
        if interval_format == "bigbed" and path.with_suffix(".bb").exists():
            return path.with_suffix(".bb")
        return path

    candidate_paths = {name: selected_interval(path) for name, path in candidate_paths.items()}
    final_paths = {name: selected_interval(path) for name, path in final_paths.items()}
    final_tss_paths = {name: selected_interval(path) for name, path in final_tss_paths.items()}
    overlap_path = selected_interval(overlap_path)
    final_state_path = selected_interval(final_state_path)
    final_tss_state_path = selected_interval(final_tss_state_path)

    assignment_path = outdir / f"{prefix}_gene_assignments.tsv"
    with assignment_path.open("w") as handle:
        handle.write(
            "gene_id\tgene_name\tchrom\tstart\tend\tintersecting_states\t"
            "candidate_sets\texcluded_from_sets\tfinal_set\tcandidate_overlap\n"
        )
        for gene in genes:
            candidates = sorted(membership.get(gene.gene_id, []))
            excluded_from = sorted(
                name for name, ids in excluded_ids_by_set.items() if gene.gene_id in ids
            )
            finals = [name for name in ordered_set_names if gene.gene_id in final_ids.get(name, set())]
            handle.write(
                f"{gene.gene_id}\t{_gene_name(gene)}\t{gene.chrom}\t{gene.start}\t{gene.end}\t"
                f"{','.join(sorted(state_by_gene[gene.gene_id]))}\t"
                f"{','.join(candidates)}\t{','.join(excluded_from)}\t"
                f"{finals[0] if finals else ''}\t"
                f"{'yes' if gene.gene_id in overlapping_ids else 'no'}\n"
            )

    pairwise_path = outdir / f"{prefix}_pairwise_overlap.tsv"
    with pairwise_path.open("w") as handle:
        handle.write("set_a\tset_b\tset_a_candidates\tset_b_candidates\tshared_genes\n")
        for index, left in enumerate(rules):
            for right in rules[index + 1 :]:
                handle.write(
                    f"{left.name}\t{right.name}\t{len(candidate_ids[left.name])}\t"
                    f"{len(candidate_ids[right.name])}\t"
                    f"{len(candidate_ids[left.name] & candidate_ids[right.name])}\n"
                )

    config_path = outdir / f"{prefix}_rules.tsv"
    with config_path.open("w") as handle:
        handle.write("set_name\tinclude_rule\texclude_if_candidate\n")
        for rule in rules:
            handle.write(
                f"{rule.name}\t{rule.expression}\t{','.join(sorted(rule.exclude_if_candidate))}\n"
            )

    summary_path = outdir / f"{prefix}_summary.tsv"
    with summary_path.open("w") as handle:
        handle.write(
            "set_name\tinclude_rule\texclude_if_candidate\tcandidate_gene_count\t"
            "overlap_removed_count\texcluded_gene_count\tfinal_gene_count\t"
            "candidate_interval\tfinal_interval\tfinal_tss_interval\n"
        )
        for rule in rules:
            handle.write(
                f"{rule.name}\t{rule.expression}\t"
                f"{','.join(sorted(rule.exclude_if_candidate))}\t"
                f"{len(candidate_ids[rule.name])}\t{len(excluded_ids_by_set[rule.name])}\t"
                f"{len(excluded_ids_by_set[rule.name])}\t{len(final_ids[rule.name])}\t"
                f"{candidate_paths[rule.name].resolve()}\t"
                f"{final_paths[rule.name].resolve()}\t"
                f"{final_tss_paths[rule.name].resolve()}\n"
            )
        if leftover_set_name:
            handle.write(
                f"{leftover_set_name}\tno candidate-set intersection\t\t{len(final_ids[leftover_set_name])}\t0\t0\t"
                f"{len(final_ids[leftover_set_name])}\t\t{final_paths[leftover_set_name].resolve()}\t"
                f"{final_tss_paths[leftover_set_name].resolve()}\n"
            )

    selected_venn_sets = list(venn_sets or [])
    if not selected_venn_sets:
        names = [rule.name for rule in rules]
        preferred_three = ["active_genes", "weak_genes", "repressed_genes"]
        preferred_two = ["active", "repressed"]
        if set(preferred_three) <= set(names):
            selected_venn_sets = preferred_three
        elif set(preferred_two) <= set(names):
            selected_venn_sets = preferred_two
        else:
            selected_venn_sets = names[: min(3, len(names))]
    if len(selected_venn_sets) not in (2, 3):
        raise ValueError("--venn-sets requires exactly two or three set names")
    unknown = [name for name in selected_venn_sets if name not in candidate_ids]
    if unknown:
        raise ValueError("Unknown Venn set names: " + ", ".join(unknown))
    venn_path = outdir / f"{prefix}_{'_vs_'.join(map(safe_name, selected_venn_sets))}_venn.png"
    venn_path = create_venn(
        selected_venn_sets,
        candidate_ids,
        venn_path,
        "Candidate gene-set overlap",
    )

    return {
        "summary": summary_path,
        "assignments": assignment_path,
        "pairwise": pairwise_path,
        "overlap_interval": overlap_path,
        "overlap_bed": overlap_path,
        "final_state_interval": final_state_path,
        "final_tss_interval": final_tss_state_path,
        "venn": venn_path,
        "rules": config_path,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite gene-sets",
        description=(
            "Create candidate gene sets from chromatin-state intersections and resolve "
            "them into mutually exclusive final categories."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("--genes-bed", required=True, help="BED3+ file containing one interval per gene.")
    parser.add_argument("--states-bed", required=True, help="BED3+ chromatin-state segmentation.")
    parser.add_argument(
        "--config",
        help=(
            "TSV with set_name and include_rule columns. An optional "
            "exclude_if_candidate column lists candidate sets whose membership excludes "
            "a gene from the current final set."
        ),
    )
    parser.add_argument(
        "--gene-set",
        action="append",
        default=[],
        metavar="NAME=RULE",
        help="Inline inclusion rule; repeat for each set. Operators: &, | and parentheses.",
    )
    parser.add_argument(
        "--leftover-set-name",
        default=None,
        help="Optional final category containing genes that belong to none of the configured candidate sets.",
    )
    parser.add_argument("--gene-id-column", type=int, default=4, help="One-based gene identifier column.")
    parser.add_argument("--state-label-column", type=int, default=4, help="One-based state label column.")
    parser.add_argument("--chrom-sizes", help="Chromosome-size table, BAM or CRAM used to filter/clip inputs and create bigBed output.")
    parser.add_argument(
        "--blacklist-bed",
        help=(
            "BED blacklist. Genes whose one-base TSS anchor overlaps the blacklist "
            "are excluded; unstranded records use their complete interval."
        ),
    )
    parser.add_argument(
        "--interval-format", choices=INTERVAL_FORMATS, default="bed",
        help="Write gene-set intervals as BED, bigBed, or both.",
    )
    parser.add_argument("--output-dir", default=".", help="Output directory (default: current directory).")
    parser.add_argument("--output-prefix", help="Output-file prefix. Default: genes and states basenames plus _gene_sets.")
    parser.add_argument(
        "--prefix-member-files",
        action="store_true",
        help=(
            "Prefix candidate/final/TSS member BED filenames with --output-prefix. "
            "The suites use this only for randomized-control filenames."
        ),
    )
    parser.add_argument(
        "--venn-sets",
        nargs="+",
        help="Exactly two or three candidate set names to display in the Venn diagram.",
    )
    add_parallel_arguments(
        parser,
        cores_option="--memory-intensive-analysis-cores",
        cores_help=(
            "Concurrent contig workers for this memory-intensive interval analysis "
            "(default: 1; independent of suite --cores)."
        ),
    )
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def _run_serial(args: argparse.Namespace) -> int:
    reporter = ProgressReporter("gene-sets")
    reporter.stage("Loading chromosome sizes, rules, and genes")
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    rules = load_rules(args.config, args.gene_set)
    genes = read_genes(args.genes_bed, args.gene_id_column, chrom_sizes)
    blacklist = load_blacklist_unbounded(args.blacklist_bed)
    genes, blacklisted_gene_anchors = filter_blacklisted_gene_anchors(
        genes, blacklist
    )
    if not genes:
        raise ValueError("No genes remained after blacklist-anchor filtering")
    reporter.stage(
        f"Loaded {len(genes):,} genes; loading chromatin-state annotations"
    )
    states = read_states(args.states_bed, args.state_label_column, chrom_sizes)
    reporter.stage(
        f"Assigning genes using {len(states):,} chromatin-state intervals"
    )
    outputs = build_gene_sets(
        genes,
        states,
        rules,
        args.output_dir,
        args.output_prefix,
        args.venn_sets,
        args.interval_format,
        chrom_sizes,
        args.leftover_set_name,
        args.prefix_member_files,
    )
    if blacklist is not None:
        metadata_path = (
            Path(args.output_dir)
            / f"{safe_name(args.output_prefix)}_blacklist_summary.tsv"
        )
        with metadata_path.open("wt", encoding="utf-8") as handle:
            handle.write("metric\tvalue\n")
            handle.write(f"blacklist_bed\t{blacklist.path}\n")
            handle.write(
                "blacklisted_gene_anchors_excluded\t"
                f"{blacklisted_gene_anchors}\n"
            )
            handle.write(f"retained_genes\t{len(genes)}\n")
        print(
            "Blacklisted gene anchors excluded: "
            f"{blacklisted_gene_anchors:,}; summary: {metadata_path}"
        )
    print(f"Gene-set summary: {outputs['summary']}")
    print(f"Final state-labelled interval: {outputs['final_state_interval']}")
    print(f"Candidate overlap Venn diagram: {outputs['venn']}")
    print(f"Genes shared by candidate sets: {outputs['overlap_interval']}")
    return 0


def run(args: argparse.Namespace) -> int:
    if not args.output_prefix:
        from nucleosuite.output_naming import input_basename
        args.output_prefix = f"{input_basename(args.genes_bed)}_{input_basename(args.states_bed)}_gene_sets"
    return run_partitioned_command(
        "gene-sets",
        args,
        _run_serial,
        runner_module="nucleosuite.gene_sets",
        runner_function="_run_serial",
        primary_attr="genes_bed",
        output_prefix_attr="output_prefix",
        output_dir_attr="output_dir",
        path_attrs=("states_bed", "blacklist_bed"),
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
