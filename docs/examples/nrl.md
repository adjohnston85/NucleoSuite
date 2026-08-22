# `nrl` example outputs

[Command documentation](../commands/nrl.md) · [Back to the README workflow](../../README.md#clickable-workflow-test-dyads--dac--nrl)

`nucleosuite nrl` detects recurring peaks in a DAC or DCC profile and fits their positions to estimate the repeat period. For dyad DAC, this period can be interpreted as nucleosome repeat length.

## Example command

```bash
nucleosuite nrl sample_145_147_dac.tsv \
  --peak-resolution 160 \
  --output-prefix sample_145_147_nrl
```

## Example outputs

`nrl` writes two principal plots:

- a profile plot showing the input profile, smoothing scales, and retained periodic peaks;
- a regression plot showing peak number against peak distance and the fitted repeat period.

A representative NRL profile/regression image can be added to this page later without changing the clickable README workflow.
