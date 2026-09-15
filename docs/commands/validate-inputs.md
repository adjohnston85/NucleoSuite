# `nucleosuite validate-inputs`

## What this command does

`validate-inputs` checks whether the main files for a NucleoSuite analysis are internally usable and compatible with one another.

## Why use it

Run it before a long multicontig or suite job so input problems are found before expensive analysis begins.

## Basic usage

```bash
nucleosuite validate-inputs \
  --bam sample.bam \
  --fasta genome.fa
```

Add the same resource/reference selections you intend to use in the final workflow when the validator supports them.

## What it checks

Depending on the supplied inputs, validation can check:

- file existence and readability;
- BAM/CRAM header and index requirements;
- FASTA/index availability;
- contig-name and contig-length compatibility;
- requested contig resolution; and
- compatibility of supporting resources used by the planned analysis.

The exact set of checks depends on which input types are supplied.

## How to use the result

A successful validation means the checked structural and reference requirements passed. Biological parameters such as fragment ranges and tissue/resource choices remain experiment-specific analysis decisions.

[Back to the command reference](../COMMAND_REFERENCE.md)
