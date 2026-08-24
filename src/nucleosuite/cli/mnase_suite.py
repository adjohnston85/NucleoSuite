"""Launch and parallelize the packaged MNase full-suite workflow."""

from __future__ import annotations

import glob
import gzip
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.resources import as_file, files
from pathlib import Path
from collections.abc import Sequence

from nucleosuite.cli.suite_routing import (
    ensure_bam_index as _ensure_bam_index,
    route_bams_by_contig as _route_bams_by_contig,
)



def _consume_values(argv: Sequence[str], option: str) -> list[str]:
    values: list[str] = []
    index = 0
    while index < len(argv):
        token = argv[index]
        if token == option:
            index += 1
            while index < len(argv) and not argv[index].startswith("-"):
                values.append(argv[index])
                index += 1
            continue
        if token.startswith(option + "="):
            values.append(token.split("=", 1)[1])
        index += 1
    return values


def _single_value(argv: Sequence[str], option: str, default: str | None = None) -> str | None:
    values = _consume_values(argv, option)
    return values[-1] if values else default


def _extract_cores(argv: Sequence[str]) -> tuple[int, list[str]]:
    output: list[str] = []
    cores = 1
    index = 0
    while index < len(argv):
        token = argv[index]
        if token == "--cores":
            if index + 1 >= len(argv):
                raise ValueError("--cores requires an integer")
            cores = int(argv[index + 1])
            index += 2
            continue
        if token.startswith("--cores="):
            cores = int(token.split("=", 1)[1])
            index += 1
            continue
        output.append(token)
        index += 1
    if cores < 1:
        raise ValueError("--cores must be at least 1")
    return cores, output


def _extract_phase_cores(
    argv: Sequence[str], option: str, default: int
) -> tuple[int, list[str]]:
    output: list[str] = []
    value = int(default)
    index = 0
    while index < len(argv):
        token = argv[index]
        if token == option:
            if index + 1 >= len(argv):
                raise ValueError(f"{option} requires an integer")
            value = int(argv[index + 1])
            index += 2
            continue
        if token.startswith(option + "="):
            value = int(token.split("=", 1)[1])
            index += 1
            continue
        output.append(token)
        index += 1
    if value < 1:
        raise ValueError(f"{option} must be at least 1")
    return value, output


def _extract_analysis_scope(argv: Sequence[str]) -> tuple[str, list[str]]:
    """Consume the multicontig analysis-scope wrapper option.

    ``combined-only`` is the default: workers create only combine prerequisites,
    and downstream analyses run once after those prerequisites are combined.
    """
    output: list[str] = []
    scope = "combined-only"
    index = 0
    while index < len(argv):
        token = argv[index]
        if token == "--analysis-scope":
            if index + 1 >= len(argv):
                raise ValueError(
                    "--analysis-scope requires combined-only or per-contig-and-combined"
                )
            scope = argv[index + 1].strip().lower()
            index += 2
            continue
        if token.startswith("--analysis-scope="):
            scope = token.split("=", 1)[1].strip().lower()
            index += 1
            continue
        output.append(token)
        index += 1
    aliases = {
        "combined": "combined-only",
        "combine-first": "combined-only",
        "all": "per-contig-and-combined",
    }
    scope = aliases.get(scope, scope)
    if scope not in {"combined-only", "per-contig-and-combined"}:
        raise ValueError(
            "--analysis-scope must be combined-only or per-contig-and-combined"
        )
    return scope, output


def _extract_combine_bigwig_method(argv: Sequence[str]) -> tuple[str, list[str]]:
    output: list[str] = []
    method = "direct"
    index = 0
    while index < len(argv):
        token = argv[index]
        if token == "--combine-bigwig-method":
            if index + 1 >= len(argv):
                raise ValueError("--combine-bigwig-method requires direct or bedgraph")
            method = argv[index + 1].strip().lower()
            index += 2
            continue
        if token.startswith("--combine-bigwig-method="):
            method = token.split("=", 1)[1].strip().lower()
            index += 1
            continue
        output.append(token)
        index += 1
    if method == "bedgraphs":
        method = "bedgraph"
    if method not in {"direct", "bedgraph"}:
        raise ValueError("--combine-bigwig-method must be direct or bedgraph")
    return method, output


def _replace_multi_option(argv: Sequence[str], option: str, values: Sequence[str]) -> list[str]:
    output: list[str] = []
    index = 0
    while index < len(argv):
        token = argv[index]
        if token == option:
            index += 1
            while index < len(argv) and not argv[index].startswith("-"):
                index += 1
            continue
        if token.startswith(option + "="):
            index += 1
            continue
        output.append(token)
        index += 1
    output.extend([option, *values])
    return output


def _replace_single_option(argv: Sequence[str], option: str, value: str) -> list[str]:
    output: list[str] = []
    index = 0
    while index < len(argv):
        token = argv[index]
        if token == option:
            index += 2
            continue
        if token.startswith(option + "="):
            index += 1
            continue
        output.append(token)
        index += 1
    output.extend([option, value])
    return output


def _remove_option(argv: Sequence[str], option: str, *, multi: bool = False) -> list[str]:
    output: list[str] = []
    index = 0
    while index < len(argv):
        token = argv[index]
        if token == option:
            index += 1
            if multi:
                while index < len(argv) and not argv[index].startswith("-"):
                    index += 1
            elif index < len(argv):
                index += 1
            continue
        if token.startswith(option + "="):
            index += 1
            continue
        output.append(token)
        index += 1
    return output


def _has_flag(argv: Sequence[str], option: str) -> bool:
    return option in argv or any(token.startswith(option + "=") for token in argv)


def _remove_flag(argv: Sequence[str], option: str) -> list[str]:
    return [token for token in argv if token != option and not token.startswith(option + "=")]


_RANDOMIZED_SUFFIX = "_randomized_control"


def _worker_sample_name(sample_name: str, contig: str) -> str:
    """Insert the contig before the randomized marker, keeping it terminal."""
    if sample_name.endswith(_RANDOMIZED_SUFFIX):
        base = sample_name[: -len(_RANDOMIZED_SUFFIX)]
        return f"{base}_{contig}{_RANDOMIZED_SUFFIX}"
    return f"{sample_name}_{contig}"


def _combine_randomized_fragments(
    worker_roots: Sequence[Path],
    combined_root: Path,
    sample_name: str,
    contigs: Sequence[str],
) -> Path:
    """Atomically combine validated per-contig randomized BED.gz files."""
    selected: list[Path] = []
    qc_inputs: list[Path] = []
    relocation_inputs: list[Path] = []
    for root in worker_roots:
        matches = sorted((root / "00_setup").glob("*.randomized.fragments.bed.gz"))
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one randomized fragment BED in {root / '00_setup'}, found {len(matches)}"
            )
        selected.append(matches[0])
        qc_matches = sorted((root / "00_setup").glob("*.randomization_qc.tsv"))
        relocation_matches = sorted(
            (root / "00_setup").glob("*.relocation_distances.tsv")
        )
        if len(qc_matches) != 1 or len(relocation_matches) != 1:
            raise RuntimeError(
                f"Expected one randomization QC table and one relocation table in "
                f"{root / '00_setup'}"
            )
        qc_inputs.append(qc_matches[0])
        relocation_inputs.append(relocation_matches[0])
    setup = combined_root / "00_setup"
    setup.mkdir(parents=True, exist_ok=True)
    output = setup / f"{sample_name}.randomized.fragments.bed.gz"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".partial", dir=setup
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    qc_output = setup / f"{sample_name}.randomization_qc.tsv"
    relocation_output = setup / f"{sample_name}.relocation_distances.tsv"
    from nucleosuite.plotting import plot_path
    relocation_plot = plot_path(setup / f"{sample_name}.relocation_distances.png")
    temporary_companions: list[tuple[Path, Path]] = []
    expected_order = {chrom: index for index, chrom in enumerate(contigs)}
    previous: tuple[int, int, int] | None = None
    count = 0
    try:
        with gzip.open(temporary, "wt", encoding="utf-8") as destination:
            for source in selected:
                with gzip.open(source, "rt", encoding="utf-8") as handle:
                    for line_number, raw in enumerate(handle, 1):
                        fields = raw.rstrip("\n").split("\t")
                        if len(fields) < 3:
                            raise RuntimeError(f"Invalid randomized BED3 line in {source}:{line_number}")
                        chrom = fields[0]
                        try:
                            start, end = int(fields[1]), int(fields[2])
                        except ValueError as exc:
                            raise RuntimeError(f"Invalid randomized coordinates in {source}:{line_number}") from exc
                        if chrom not in expected_order or start < 0 or end <= start:
                            raise RuntimeError(f"Invalid randomized interval in {source}:{line_number}")
                        key = (expected_order[chrom], start, end)
                        if previous is not None and key < previous:
                            raise RuntimeError("Per-contig randomized fragments are not globally sorted")
                        previous = key
                        destination.write(f"{chrom}\t{start}\t{end}\n")
                        count += 1
        if count == 0:
            raise RuntimeError("Combined randomized fragment BED would be empty")
        with gzip.open(temporary, "rt", encoding="utf-8") as check:
            validated = sum(1 for _ in check)
        if validated != count:
            raise RuntimeError("Combined randomized fragment BED failed record-count validation")
        from nucleosuite.combine import _combine_generic_tsv, _combine_randomization_qc
        from nucleosuite.profile_plots import plot_count_profile

        for published, suffix in (
            (qc_output, ".tsv"),
            (relocation_output, ".tsv"),
            (relocation_plot, relocation_plot.suffix),
        ):
            descriptor, name = tempfile.mkstemp(
                prefix=f".{published.stem}.", suffix=suffix, dir=setup
            )
            os.close(descriptor)
            temporary_companions.append((Path(name), published))
        temporary_qc, temporary_relocation, temporary_plot = (
            item[0] for item in temporary_companions
        )
        _combine_randomization_qc(qc_inputs, temporary_qc)
        _combine_generic_tsv(relocation_inputs, temporary_relocation)
        temporary_plot = plot_count_profile(
            temporary_relocation,
            temporary_plot,
            x_column="relocation_bp",
            y_column="count",
            xlabel="Fragment relocation (bp)",
            ylabel="Fragment count",
            title="Randomized fragment relocation distances",
            vertical_zero=True,
        )
        if any(
            not path.is_file() or path.stat().st_size == 0
            for path, _ in temporary_companions
        ):
            raise RuntimeError("Combined randomization QC output validation failed")
        for path, published in temporary_companions:
            os.replace(path, published)
        # Publish the validated fragment input last, after all companion QC files.
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
        for path, _published in temporary_companions:
            path.unlink(missing_ok=True)
    return output


def _expand_bam_inputs(argv: Sequence[str]) -> list[str]:
    """Expand BAM paths and globs once before contig workers are launched."""

    expanded: list[str] = []
    for token in _consume_values(argv, "--bam"):
        matches = [token] if Path(token).is_file() else glob.glob(token)
        if not matches:
            raise ValueError(f"BAM input did not match any files: {token}")
        for match in matches:
            path = Path(match).resolve()
            if path.suffix.lower() != ".bam":
                raise ValueError(f"BAM input is not a .bam file: {path}")
            expanded.append(str(path))
    return sorted(dict.fromkeys(expanded), key=_natural_path_key)


def _expand_fragment_inputs(argv: Sequence[str]) -> list[str]:
    """Expand fragment paths and globs once for stable provenance and routing."""
    expanded: list[str] = []
    for token in _consume_values(argv, "--fragments"):
        matches = [token] if Path(token).is_file() else glob.glob(token)
        if not matches:
            raise ValueError(f"Fragment input did not match any files: {token}")
        expanded.extend(
            str(Path(match).resolve())
            for match in matches
            if Path(match).is_file()
        )
    if not expanded:
        raise ValueError("No fragment inputs were found")
    return sorted(dict.fromkeys(expanded), key=_natural_path_key)


def _natural_path_key(value: str) -> list[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value)]


def _derive_sample_name(
    input_paths: Sequence[str],
    explicit: str | None = None,
    *,
    fallback: str = "multi_bam",
) -> str:
    """Derive the chromosome-independent sample prefix used by all workers."""

    if explicit:
        cleaned = re.sub(r"[^A-Za-z0-9._-]", "_", explicit)
        if not cleaned:
            raise ValueError("--sample-name does not contain any usable characters")
        return cleaned
    stems = []
    for value in input_paths:
        name = Path(value).name
        for suffix in (".bed.gz", ".bigBed", ".bigbed", ".bam", ".bed", ".bb"):
            if name.lower().endswith(suffix.lower()):
                name = name[: -len(suffix)]
                break
        stems.append(name)
    if len(stems) == 1:
        return re.sub(r"[^A-Za-z0-9._-]", "_", stems[0])
    stripped = [
        re.sub(r"[._-]chr(?:\d+|X|Y|M|MT)$", "", stem, flags=re.IGNORECASE)
        for stem in stems
    ]
    candidate = (
        stripped[0]
        if stripped and all(value == stripped[0] for value in stripped)
        else fallback
    )
    return re.sub(r"[^A-Za-z0-9._-]", "_", candidate)


def _resolve_contigs(argv: Sequence[str]) -> tuple[list[str], list[tuple[str, int]], list[tuple[str, int]]]:
    """Resolve requested contigs in the BAM-derived output namespace."""
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("mnase-suite multicontig processing requires pysam") from exc

    bam_paths = _expand_bam_inputs(argv) if _consume_values(argv, "--bam") else []
    if bam_paths:
        from nucleosuite.core.bam_headers import merge_bam_reference_headers_with_aliases
        handles = [pysam.AlignmentFile(path, "rb") for path in bam_paths]
        try:
            merged = merge_bam_reference_headers_with_aliases(handles)
            all_rows = list(zip(merged.references, merged.lengths))
        finally:
            for handle in handles:
                handle.close()
    else:
        fasta_path = _single_value(argv, "--fasta")
        if not fasta_path:
            raise ValueError("--fasta is required")
        try:
            fasta = pysam.FastaFile(fasta_path)
        except (OSError, ValueError):
            pysam.faidx(fasta_path)
            fasta = pysam.FastaFile(fasta_path)
        try:
            all_rows = [
                (chrom, int(fasta.get_reference_length(chrom)))
                for chrom in fasta.references
            ]
        finally:
            fasta.close()

    from nucleosuite.core.regions import expand_contig_tokens
    references = [chrom for chrom, _length in all_rows]
    lengths = dict(all_rows)
    tokens = _consume_values(argv, "--contigs") or ["autosomes"]
    specs = expand_contig_tokens(tokens, references)
    if any(":" in spec for spec in specs):
        raise ValueError("mnase-suite multicore mode accepts whole contigs, not coordinate ranges")
    selected: list[str] = []
    for chrom in specs:
        if chrom not in selected:
            selected.append(chrom)
    return selected, [(chrom, lengths[chrom]) for chrom in selected], all_rows


def _run_shell(
    script_path: str,
    argv: list[str],
    log_path: str | None = None,
    env_updates: dict[str, str] | None = None,
) -> int:
    environment = os.environ.copy()
    if env_updates:
        environment.update(env_updates)
    if log_path is None:
        return int(
            subprocess.run(
                ["bash", script_path, *argv], check=False, env=environment
            ).returncode
        )
    path = Path(log_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            ["bash", script_path, *argv],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            env=environment,
        )
    return int(completed.returncode)


def _parallel_main(
    script_path: str,
    argv: list[str],
    cores: int,
    streaming_combine_cores: int,
    indexed_combine_cores: int,
    combine_chunk_bp: int,
    analysis_cores: int,
    memory_intensive_analysis_cores: int,
    combine_bigwig_method: str,
    analysis_scope: str,
) -> int:
    randomized_mode = _has_flag(argv, "--randomize")
    contigs, chrom_sizes, all_chrom_sizes = _resolve_contigs(argv)
    bam_paths = _expand_bam_inputs(argv) if _consume_values(argv, "--bam") else []
    explicit_sample = _single_value(argv, "--sample-name")
    fragment_paths = (
        _expand_fragment_inputs(argv) if _consume_values(argv, "--fragments") else []
    )
    global_sample = (
        _derive_sample_name(bam_paths, explicit_sample)
        if bam_paths
        else _derive_sample_name(
            fragment_paths, explicit_sample, fallback="multi_fragments"
        )
        if fragment_paths
        else explicit_sample
    )
    if randomized_mode and global_sample and not global_sample.endswith(_RANDOMIZED_SUFFIX):
        global_sample += _RANDOMIZED_SUFFIX
    if bam_paths:
        argv = _replace_multi_option(argv, "--bam", bam_paths)
    if fragment_paths:
        argv = _replace_multi_option(argv, "--fragments", fragment_paths)
    if global_sample:
        argv = _replace_single_option(argv, "--sample-name", global_sample)
    if len(contigs) <= 1:
        return _run_shell(
            script_path,
            argv,
            env_updates={
                "NUCLEOSUITE_SUITE_CORES": str(analysis_cores),
                "NUCLEOSUITE_MEMORY_ANALYSIS_CORES": str(
                    memory_intensive_analysis_cores
                ),
            },
        )
    outdir_value = _single_value(argv, "--outdir")
    if not outdir_value:
        raise ValueError("--outdir is required")
    root = Path(outdir_value).resolve()
    per_root = root / "per_contig"
    combined_root = root / "combined"
    per_root.mkdir(parents=True, exist_ok=True)
    combined_root.mkdir(parents=True, exist_ok=True)

    from nucleosuite.core.chrom_sizes import write_chrom_sizes_table
    suite_sample = global_sample or "sample"
    support_prefix = f"{suite_sample}_" if randomized_mode else ""
    shared_setup = root / "00_setup"
    shared_sizes = write_chrom_sizes_table(
        all_chrom_sizes, shared_setup / f"{support_prefix}analysis.chrom.sizes"
    )
    argv = _replace_single_option(
        argv, "--analysis-chrom-sizes-source", str(shared_sizes)
    )

    routed_bams = _route_bams_by_contig(bam_paths, contigs) if bam_paths else {contig: [] for contig in contigs}
    active_contigs = [contig for contig in contigs if not bam_paths or routed_bams[contig]]
    skipped_contigs = [contig for contig in contigs if bam_paths and not routed_bams[contig]]
    for contig in skipped_contigs:
        print(f"Skipping {contig}: no mapped reads were found in the supplied BAM files")
    if not active_contigs:
        raise ValueError("No mapped reads were found for any requested contig")
    active_sizes = [(contig, length) for contig, length in chrom_sizes if contig in set(active_contigs)]

    failures: list[str] = []
    with ProcessPoolExecutor(max_workers=min(cores, len(active_contigs))) as executor:
        futures = {}
        for contig in active_contigs:
            safe = contig.replace("/", "_").replace("\\", "_").replace(":", "_")
            contig_out = per_root / safe
            worker_args = _replace_multi_option(argv, "--contigs", [contig])
            worker_args = _replace_single_option(worker_args, "--outdir", str(contig_out))
            if bam_paths:
                worker_args = _replace_multi_option(worker_args, "--bam", routed_bams[contig])
            worker_sample = _worker_sample_name(suite_sample, contig)
            if global_sample:
                worker_args = _replace_single_option(worker_args, "--sample-name", worker_sample)
            # In the default combine-first mode, workers create only the
            # active-mode tracks and sufficient-statistic outputs needed
            # for combination. All analytical stages run once after combination.
            if analysis_scope == "combined-only":
                worker_args.append("--combine-prerequisites-only")
            else:
                # Expression analyses use the combined output.
                if _single_value(argv, "--expression"):
                    worker_args.append("--skip-gene-expression")
                if "--skip-tss-expression-quintiles" not in worker_args:
                    worker_args.append("--skip-tss-expression-quintiles")
            worker_environment = {
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "NUMEXPR_NUM_THREADS": "1",
            }
            if combine_bigwig_method == "bedgraph":
                worker_environment.update({
                    "NUCLEOSUITE_STAGED_BEDGRAPH_ROOT": str(
                        combined_root / "temporary_bedgraph_combine"
                    ),
                    "NUCLEOSUITE_STAGED_BEDGRAPH_SOURCE_ROOT": str(contig_out),
                    "NUCLEOSUITE_STAGED_BEDGRAPH_SOURCE_ID": safe,
                })
            future = executor.submit(
                _run_shell,
                script_path,
                worker_args,
                str(
                    contig_out
                    / (
                        f"{worker_sample}_mnase-suite.log"
                        if randomized_mode
                        else "mnase-suite.log"
                    )
                ),
                worker_environment,
            )
            futures[future] = contig
        for future in as_completed(futures):
            contig = futures[future]
            try:
                code = future.result()
            except Exception as error:
                failures.append(f"{contig}: {error}")
                continue
            if code:
                failures.append(f"{contig}: exit code {code}")
            else:
                print(f"Completed MNase suite for {contig}")
    if failures:
        print(
            "WARNING: Some per-contig MNase workflows recorded failures. "
            "Available upstream outputs will still be combined."
        )
        for failure in failures:
            print(f"WARNING: {failure}")

        if analysis_scope == "combined-only":
            raise RuntimeError(
                "At least one per-contig prerequisite workflow failed; no partial "
                "combined analysis was created: " + "; ".join(failures)
            )

    from nucleosuite.combine import combine_directory_trees

    worker_roots = [
        per_root / contig.replace("/", "_").replace("\\", "_").replace(":", "_")
        for contig in active_contigs
    ]
    upstream_roots = ("01_combined_tracks",)
    result = combine_directory_trees(
        worker_roots,
        combined_root,
        chrom_sizes=active_sizes,
        include_roots=upstream_roots,
        exclude_parts=("compare_positions", "score_frequencies"),
        sample_name=global_sample,
        bigwig_method=combine_bigwig_method,
        strict_complete=True,
        cores=streaming_combine_cores,
        streaming_cores=streaming_combine_cores,
        indexed_cores=indexed_combine_cores,
        bigwig_chunk_size=combine_chunk_bp,
    )
    setup_dir = combined_root / "00_setup"
    setup_dir.mkdir(parents=True, exist_ok=True)
    combined_contigs_table = setup_dir / f"{support_prefix}combined_chromosomes.tsv"
    with combined_contigs_table.open("w", encoding="utf-8") as handle:
        handle.write("chromosome\tlength\tstatus\n")
        for contig, length in active_sizes:
            handle.write(f"{contig}\t{length}\tincluded\n")
        for contig in skipped_contigs:
            length = dict(chrom_sizes).get(contig, "")
            handle.write(f"{contig}\t{length}\tskipped_no_mapped_reads\n")

    for warning in result["warnings"]:
        print(f"WARNING: {warning}")
    combine_log = Path(
        result.get("combine_log", combined_root / f"{support_prefix}combine_steps.log")
    )
    print(f"Combine log: {combine_log}")
    bigwig_warnings = [warning for warning in result["warnings"] if ".bw" in warning or ".bigWig" in warning]
    if bigwig_warnings:
        raise RuntimeError(
            "Combined BigWig creation failed. Tabular outputs were retained. "
            f"See {combine_log} and rerun the workflow after resolving the reported error."
        )

    randomized_bed: Path | None = None
    if randomized_mode:
        if not global_sample:
            raise RuntimeError("Randomized multicontig runs require a resolvable sample name")
        randomized_bed = _combine_randomized_fragments(
            worker_roots, combined_root, global_sample, active_contigs
        )

    # The final serial pass reuses the combined tracks and interval calls, and
    # runs analytical stages once across exactly the chromosomes that
    # contributed complete prerequisite outputs.
    final_args = _replace_single_option(argv, "--outdir", str(combined_root))
    final_args = _replace_multi_option(final_args, "--contigs", active_contigs)
    if randomized_bed is not None:
        final_args = _remove_option(final_args, "--bam", multi=True)
        final_args = _remove_option(final_args, "--fragments", multi=True)
        final_args = _remove_flag(final_args, "--randomize")
        final_args.extend(["--fragments", str(randomized_bed), "--randomized-control-input"])
        for path in bam_paths:
            final_args.extend(["--provenance-bam", path])
        for path in fragment_paths:
            final_args.extend(["--provenance-fragment", path])
    final_args.append("--resume")
    final_args.append("--trusted-combined-prerequisites")
    combined_log = combined_root / (
        f"{suite_sample}_mnase-suite-combined.log"
        if randomized_mode
        else "mnase-suite-combined.log"
    )
    code = _run_shell(
        script_path,
        final_args,
        str(combined_log),
        {
            "NUCLEOSUITE_SUITE_CORES": str(analysis_cores),
            "NUCLEOSUITE_MEMORY_ANALYSIS_CORES": str(
                memory_intensive_analysis_cores
            ),
        },
    )
    if code:
        failures.append(f"combined workflow: exit code {code}; see {combined_log}")
    print(f"Per-contig outputs: {per_root}")
    print(f"Combined outputs: {combined_root}")
    if failures:
        raise RuntimeError(
            "The MNase suite attempted all remaining steps but recorded failures: "
            + "; ".join(failures)
        )
    return 0


def validate_argv(argv: Sequence[str] | None = None) -> None:
    """Validate suite arguments without creating outputs or starting a job."""
    args = list(argv or [])
    from nucleosuite.cli.suite_paired import extract_paired_options
    paired, _fdr, args = extract_paired_options(args)
    if paired and _has_flag(args, "--randomize"):
        raise ValueError("--with-randomized-control cannot be combined with --randomize")
    from nucleosuite.plotting import extract_plotting_argv
    args, plot_env = extract_plotting_argv(args)
    if plot_env:
        os.environ.update(plot_env)
    cores, without_cores = _extract_cores(args)
    combine_cores, without_combine_cores = _extract_phase_cores(
        without_cores, "--combine-cores", cores
    )
    streaming_combine_cores, without_streaming_combine = _extract_phase_cores(
        without_combine_cores, "--streaming-combine-cores", combine_cores
    )
    indexed_combine_cores, without_indexed_combine = _extract_phase_cores(
        without_streaming_combine, "--indexed-combine-cores", 1
    )
    combine_chunk_bp, without_combine_chunk = _extract_phase_cores(
        without_indexed_combine, "--combine-chunk-bp", 100_000
    )
    analysis_cores, without_analysis_cores = _extract_phase_cores(
        without_combine_chunk, "--analysis-cores", cores
    )
    memory_analysis_cores, without_memory_analysis = _extract_phase_cores(
        without_analysis_cores, "--memory-intensive-analysis-cores", 1
    )
    if streaming_combine_cores > cores or analysis_cores > cores:
        raise ValueError(
            "--streaming-combine-cores and --analysis-cores cannot exceed --cores"
        )
    _scope, without_scope = _extract_analysis_scope(without_memory_analysis)
    _method, shell_args = _extract_combine_bigwig_method(without_scope)
    script_resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    with as_file(script_resource) as script_path:
        completed = subprocess.run(
            ["bash", str(script_path), *shell_args, "--validate-only"],
            capture_output=True,
            text=True,
            check=False,
        )
    if completed.returncode:
        if completed.stdout:
            sys.stdout.write(completed.stdout)
        if completed.stderr:
            sys.stderr.write(completed.stderr)
        raise SystemExit(int(completed.returncode))


def main(argv: Sequence[str] | None = None) -> int:
    """Run the bundled workflow, optionally with one worker per contig."""
    args = list(argv or [])
    from nucleosuite.cli.suite_paired import (
        annotate_suite_combined_peaks,
        extract_paired_options,
    )
    paired, paired_fdr, args = extract_paired_options(args)
    if paired:
        if _has_flag(args, "--randomize"):
            raise ValueError("--with-randomized-control cannot be combined with --randomize")
        outdir = _single_value(args, "--outdir")
        if not outdir:
            raise ValueError("--with-randomized-control requires --outdir")
        args = _replace_single_option(args, "--interval-format", "both")
        observed_code = main(args)
        if observed_code:
            return int(observed_code)
        randomized_code = main([*args, "--randomize"])
        if randomized_code:
            return int(randomized_code)
        if "--dry-run" in args:
            return 0
        outputs = annotate_suite_combined_peaks(
            outdir,
            suite_name="mnase-suite",
            fdr_threshold=paired_fdr,
        )
        for label, result in outputs.items():
            print(f"{label}_peaks_empirical_fdr\t{result.annotated_path}")
            if result.significant_path is not None:
                print(f"{label}_peaks_significant\t{result.significant_path}")
        return 0
    from nucleosuite.plotting import extract_plotting_argv
    args, plot_env = extract_plotting_argv(args)
    if plot_env:
        os.environ.update(plot_env)
    cores, without_cores = _extract_cores(args)
    combine_cores, without_combine_cores = _extract_phase_cores(
        without_cores, "--combine-cores", cores
    )
    streaming_combine_cores, without_streaming_combine = _extract_phase_cores(
        without_combine_cores, "--streaming-combine-cores", combine_cores
    )
    indexed_combine_cores, without_indexed_combine = _extract_phase_cores(
        without_streaming_combine, "--indexed-combine-cores", 1
    )
    combine_chunk_bp, without_combine_chunk = _extract_phase_cores(
        without_indexed_combine, "--combine-chunk-bp", 100_000
    )
    analysis_cores, without_analysis_cores = _extract_phase_cores(
        without_combine_chunk, "--analysis-cores", cores
    )
    memory_analysis_cores, without_memory_analysis = _extract_phase_cores(
        without_analysis_cores, "--memory-intensive-analysis-cores", 1
    )
    if streaming_combine_cores > cores or analysis_cores > cores:
        raise ValueError(
            "--streaming-combine-cores and --analysis-cores cannot exceed --cores"
        )
    analysis_scope, without_scope = _extract_analysis_scope(without_memory_analysis)
    combine_bigwig_method, shell_args = _extract_combine_bigwig_method(without_scope)
    script_resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    with as_file(script_resource) as script_path:
        if cores == 1 or "--help" in shell_args or "-h" in shell_args or "--help-plotting" in shell_args:
            return _run_shell(
                str(script_path), shell_args,
                env_updates={
                    "NUCLEOSUITE_SUITE_CORES": str(analysis_cores),
                    "NUCLEOSUITE_MEMORY_ANALYSIS_CORES": str(
                        memory_analysis_cores
                    ),
                },
            )
        return _parallel_main(
            str(script_path), shell_args, cores, streaming_combine_cores,
            indexed_combine_cores, combine_chunk_bp, analysis_cores,
            memory_analysis_cores, combine_bigwig_method, analysis_scope
        )
