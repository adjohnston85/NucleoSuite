# `dyads` example outputs

[Command documentation](../commands/dyads.md) · [Back to the README workflow](../../README.md#clickable-workflow-test-dyads--dac--nrl)

`nucleosuite dyads` generates a fragment-centre signal track, typically as a BigWig. That track can then be analysed with `dac`.

## Example command

```bash
nucleosuite dyads \
  --bam sample.bam \
  --frag-lower 145 \
  --frag-upper 147 \
  --output-format bigwig \
  --out-prefix sample_145_147
```

## Example signal

The left panel below illustrates regularly spaced dyad positions that produce a periodic DAC signal. A representative browser-style dyad-track image can be added here later without changing the README link.

![Periodic dyad signal and corresponding DAC](../images/dac_periodicity_example.png)

## Next step

Use the dyad BigWig as input to [`dac`](dac.md).
