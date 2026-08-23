# `nucleosuite randomize-fragments`

## What this command does

`randomize-fragments` creates a control fragment set by moving accepted fragments to new genomic coordinates while retaining chromosome and fragment length.

## Why use it

Use randomized fragments as a positional null model. Process observed and randomized fragments with the same downstream commands and settings.

The comparison measures the effect of genomic placement while retaining the selected fragment properties and placement constraints.

## Dinucleotide-matched randomization

With `--method dinucleotide`, NucleoSuite chooses one fragment boundary, records its dinucleotide, and searches a local reference block for valid alternative placements that preserve that selected boundary dinucleotide.

A valid placement must:

- preserve chromosome and fragment length;
- fit completely inside the search block;
- contain only canonical A/C/G/T sequence;
- differ from the original coordinate;
- respect the randomized-coordinate multiplicity cap; and
- avoid the effective blacklist.

If the selected boundary has no valid placement, the opposite boundary is tried before the configured fallback.

## Uniform randomization

`--method uniform` chooses another valid coordinate within the same fixed local search block without requiring a matching boundary dinucleotide. The placement still obeys the same fragment-length, sequence, original-coordinate, multiplicity, and blacklist constraints.

## Fallback behaviour

For failed dinucleotide matching:

```text
--fallback uniform   try a valid uniform placement
--fallback skip      omit the fragment
```

## Typical use

```bash
nucleosuite randomize-fragments \
  --bam sample.bam \
  --fasta genome.fa \
  --contigs chr1-22,chrX \
  --cores 8 \
  --method dinucleotide \
  --fallback uniform \
  --seed 12345 \
  --output-prefix sample_randomized
```

A fixed seed makes the randomization reproducible.

## Process the control like the observed data

For example:

```bash
nucleosuite pns \
  --fragments sample_randomized.randomized.fragments.bed.gz \
  --fasta genome.fa \
  --out-prefix sample_randomized_pns
```

The suite commands provide randomized-only execution with `--randomize`.

## What it writes

Outputs include:

- the randomized fragment BED/BED.gz;
- QC counts for matched, uniform, fallback, and skipped fragments;
- collision/multiplicity information;
- relocation-distance counts; and
- a relocation-distance figure.

## Randomization controls

- `--search-window` sets the maximum local reference block length;
- `--anchor-prob-start` controls the probability of initially matching the start boundary;
- `--max-randomized-per-coordinate` caps repeated randomized coordinates; `0` disables the cap;
- `--blacklist-bed` excludes source fragments and candidate placements overlapping listed regions.

## Plot customization

Relocation-distance figures use the shared options in [Plot customization](../PLOTTING.md).

## Automatic output naming

If `--output-prefix` is omitted, NucleoSuite derives the randomized-output prefix from the primary input basename in the current directory. An explicit `--output-prefix` overrides this default.

[Back to the command reference](../COMMAND_REFERENCE.md)
