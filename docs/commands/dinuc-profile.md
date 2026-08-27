# `nucleosuite dinuc-profile`

## What this command does

`dinuc-profile` aligns fragments by their centre and calculates how often each dinucleotide occurs at each relative position.

## Why use it

Use it to examine sequence periodicity around nucleosome-sized fragments, compare WW/SS patterns between fragment-length classes, or generate sequence profiles for cfDNA/MNase comparisons.

## How it works

Each fragment is aligned to the right-hand central base. Dinucleotides are counted by their position relative to that centre, then converted to frequencies using the number of valid contributing fragments at each position.

A fragment is included in sequence profiling only when its complete extracted sequence contains canonical A/C/G/T bases. The exact coordinate and frequency definitions are in [Dinucleotide profiles and WW/SS classes](../ALGORITHMS.md#dinucleotide-profiles-and-wwss-classes).

## Basic usage

```bash
nucleosuite dinuc-profile \
  --bam sample.bam \
  --fasta genome.fa \
  --frag-lower 145 \
  --frag-upper 145 \
  --out-prefix sample_145
```

Repeat with another exact or narrow fragment-length class when you want to compare their rotational sequence profiles.

## Outputs

The command writes per-position counts and frequencies, optional WW/SS summary profiles, and profile figures. Multicontig runs sum raw counts and valid-fragment denominators before recalculating frequencies.

## Plot customization

Profile figures use the shared plotting interface in [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
