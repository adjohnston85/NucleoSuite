# NucleoSuite 0.10.0 release notes

Version 0.10.0 is a full workflow release built from the 0.9.6 command surface. Existing command names remain available, with the additions and suite fixes below.

## New `chip-suite`

- Accepts one or more target BAMs and one or more matched-control BAMs.
- Defaults to TNS with accepted fragment lengths 120–500 bp.
- Supports BNS and PNS through `--scoring-method`.
- Defaults to bootstrap-stabilized automatic fragment-mode estimation from seeded random genomic-block sampling.
- Estimates target and control separately and uses an equal-histogram-weight pooled mode by default.
- Accepts `--mode INT` to bypass automatic estimation and use the supplied mode exactly for both samples.
- Divides each centred track by the finite, non-zero mean of its own matching positive track: `posTNS`, `posBNS`, or `posPNS`.
- Calls target and control peaks independently. The control is not subtracted from the target.
- Uses conservative target/control peak competition for empirical peak FDR and applies identical clustering rules to target and control winners for cluster FDR.

## Empirical peak FDR

The new `pns-peak-fdr` command compares observed PNS peaks with one or more identically processed fragment-randomized peak callsets. It preserves all observed BED fields and appends `empirical_fdr` as the final field.

The complete annotated BED is always written. `--fdr 0.05` adds a significant-peak BED without replacing the complete file.

## Paired observed/randomized suites

`cfdna-suite` and `mnase-suite` now accept `--with-randomized-control`. One invocation runs the complete observed workflow and the complete randomized workflow with identical settings, then annotates the observed combined nucleosome and breakpoint BEDs with empirical FDR. The existing `--randomize` option remains available for randomized-only runs.

## Suite reliability fixes carried into 0.10.0

- Both cfDNA and MNase suite scripts now pass `--output-dir "$COMBINED_TRACK_DIR"` to `nucleosuite tracks`. This keeps each manifest output prefix inside its declared `01_combined_tracks` tree during per-contig and multicore suite execution.
- The two MNase peak-output guards now call the defined `fatal` function instead of the undefined `die` function.
- Regression coverage includes static shell-contract tests and a small dependency-gated end-to-end run through `01_combined_tracks` for both suites.

## Upgrade check

After installing the 0.10.0 wheel, verify the command surface with:

```bash
nucleosuite --version
nucleosuite chip-suite --help
nucleosuite pns-peak-fdr --help
nucleosuite mnase-suite --help-all
nucleosuite cfdna-suite --help-all
```

For the full command inventory, see the [command reference](COMMAND_REFERENCE.md).
