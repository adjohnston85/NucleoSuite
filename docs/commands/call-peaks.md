# `nucleosuite call-peaks`

## What this command does

`call-peaks` applies the PNS positive-region caller or the WPS caller to an existing compatible BigWig signal. Use it when calling parameters should be changed, a previously combined signal should be called, or signal generation and peak calling should be separate stages.

## Why use it

Use it to separate signal generation from feature calling or to apply a consistent caller to an already combined track.

## Choose the caller

```text
--peak-caller pns
--peak-caller wps
```

Match `pns` to a signed PNS BigWig and `wps` to the WPS-family signal being evaluated. See [PNS peak calling](../ALGORITHMS.md#pns-peak-calling) and [WPS peak calling](../ALGORITHMS.md#wps-peak-calling) for the segmentation rules.

## PNS example

```bash
nucleosuite call-peaks \
  --input-bigwig sample_pns.bw \
  --peak-caller pns \
  --out-prefix sample_pns_calls
```

The PNS caller segments positive score regions directly. Breakpoint calls apply the same region logic to the sign-inverted signal. Text BED scores remain six-decimal floating-point values after `--peak-score-scale`.

For bigBed output, `--bigbed-score-scale` converts the BED score to the integer 0–1000 field. The PNS default is `1`, so native values are rounded and clamped without additional rescaling. An explicit multiplier overrides the default.

## WPS example

```bash
nucleosuite call-peaks \
  --input-bigwig sample_sm_mWPS.bw \
  --peak-caller wps \
  --out-prefix sample_wps_calls
```

The WPS caller expects the WPS-family signal whose positive regions and above-median subruns should be evaluated. The standalone `wps` workflow calls from `sm_mWPS` by default.

## Outputs and coordinate handling

Depending on the selected interval format, the command writes nucleosome-region and/or breakpoint-peak BED/bigBed files plus summaries and metadata describing the calling parameters. BED8 output stores the representative call centre in column 7, which can be used directly by `distances --position-column 7`.

With multicontig processing, combined interval files contain calls from all selected contigs. `--blacklist-bed` excludes called intervals that overlap selected blacklist regions.

[Back to the command reference](../COMMAND_REFERENCE.md)
