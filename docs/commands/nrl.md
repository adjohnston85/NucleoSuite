# `nucleosuite nrl`

## What this command does

`nrl` estimates a recurring period from repeated peaks in a DAC or DCC distance profile. For nucleosome dyad DAC, that recurring period can be interpreted as a nucleosome repeat length (NRL).

## Why use it

Use `nrl` to reduce a repeated DAC or DCC peak series to one period estimate that can be compared across signals or conditions.

## How it works

The long-range peak caller is controlled by one setting: `--peak-resolution`.

With the default resolution of 160 bp:

1. candidate repeat peaks must be at least 160 bp apart;
2. the profile is smoothed over 51 bp to find the broad peak locations;
3. each detected peak is refined to the strongest local maximum in a 21 bp-smoothed profile; and
4. the refined peak distances are fitted against peak number to estimate the repeating period.

The two smoothing windows are derived from the selected resolution:

- detection window = resolution / 3;
- local-maximum window = resolution / 6.

Each value is rounded **down** to the permitted `10n + 1` series: 11, 21, 31, 41, 51 bp, and so on. A derived value below 11 bp means no smoothing.

For example, with `--peak-resolution 160`:

```text
160 / 3 = 53.3  -> 51 bp detection smoothing
160 / 6 = 26.7  -> 21 bp local-maximum smoothing
```

See [Nucleosome repeat length](../ALGORITHMS.md#nucleosome-repeat-length) for the exact calculations.

## Typical use

```bash
nucleosuite nrl sample_dac.tsv \
  --peak-resolution 160 \
  --output-prefix sample_nrl
```

The default is already 160 bp, so the option can be omitted when that resolution is appropriate.

To use a different long-range resolution, change one value:

```bash
nucleosuite nrl sample_dac.tsv \
  --peak-resolution 200 \
  --output-prefix sample_nrl_resolution200
```

A 200 bp resolution gives a 61 bp detection window and a 31 bp local-maximum window.

For analyses where no resolution-based smoothing or peak separation is wanted, use:

```bash
nucleosuite nrl sample_dac.tsv \
  --peak-resolution 0 \
  --output-prefix sample_periodicity
```

The cfDNA and MNase suites use resolution 0 for their separate short-range periodicity summaries while using the configured long-range NRL resolution for the main NRL analysis.

## How to interpret the result

Suppose the retained DAC maxima are near:

```text
370, 555, 740, 925, 1110 bp
```

Their adjacent spacings are close to 185 bp. Fitting peak distance against peak number gives a slope near 185 bp per peak, which is the estimated repeating period.

`R²` reports how closely the retained peak positions follow one regularly spaced series. A high R² means the peak positions are well described by one repeat length.

## What it writes

The command writes:

- `_profile.tsv`, containing the unsmoothed, local-max-smoothed, and detection-smoothed profiles and the called-peak flags;
- `_peaks.tsv`, containing each refined peak and the broad detection peak from which it was derived;
- `_regression.tsv`, containing the selected resolution, derived smoothing windows, fitted slope, intercept, R², slope standard error, and mean adjacent peak spacing;
- `_profile.png`, showing the input profile, the two smoothing scales, and the retained peaks; and
- `_regression.png`, showing peak number against peak distance and the fitted repeat period.

## Plot customization

Retained periodic maxima are labelled by default. Labels show the peak distance or lag and are centred directly above the called peak. Use `--plot-label-points none` to hide them. See [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
