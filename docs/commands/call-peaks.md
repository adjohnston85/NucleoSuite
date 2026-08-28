# `nucleosuite call-peaks`

## What this command does

`call-peaks` applies the PNS positive-region caller or the WPS caller to an existing compatible BigWig signal.

## Why use it

Use it to change calling parameters, call a previously combined signal, or run track generation and peak calling as separate workflow stages.

## Choose the caller that matches the signal

```text
--peak-caller pns
--peak-caller wps
```

The selected caller determines how the input signal is segmented and scored. Match `pns` to an PNS signal and `wps` to the WPS-family signal being evaluated.

See [PNS peak calling](../ALGORITHMS.md#pns-peak-calling) and [WPS peak calling](../ALGORITHMS.md#wps-peak-calling) for the exact definitions.

## Basic usage

```bash
nucleosuite call-peaks \
  --input-bigwig sample_pns.bw \
  --peak-caller pns \
  --out-prefix sample_pns_calls
```

The shared nucleosome-score caller segments positive score regions directly. Breakpoint calls apply the same region logic to the sign-inverted signal.

Text BED scores remain six-decimal floats. `--bigbed-score-scale` defaults to 1; the bigBed score is rounded and clamped to its integer range of 0–1000. An explicit multiplier affects bigBed scores only.

## WPS example

```bash
nucleosuite call-peaks \
  --input-bigwig sample_sm_mWPS.bw \
  --peak-caller wps \
  --out-prefix sample_wps_calls
```

The WPS caller expects the WPS-family signal whose positive regions and above-median subruns should be evaluated. The default standalone `wps` workflow calls from `sm_mWPS`.

## Outputs

Depending on the selected interval format, the command writes nucleosome-region and/or breakpoint-peak BED/bigBed files plus summaries/metadata describing the calling parameters.

NucleoSuite BED8 interval output stores the representative call centre in column 7. That position can be used directly by `distances --position-column 7`.

## Chromosome-wise processing

With multicontig processing, the combined interval files contain the calls from all selected contigs.

## Blacklist handling

`--blacklist-bed` excludes called intervals that overlap selected blacklist regions.

[Back to the command reference](../COMMAND_REFERENCE.md)
