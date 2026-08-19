# NucleoSuite 0.8.8

- Added generic `flank-spacing` analysis for category-wise distributions of distances between nucleosomes flanking reference sites.
- Density distributions are default; raw counts are optional. Categories are ranked by the ratio of curve heights at configurable positions (190/260 bp by default), with the lowest ratio ranked first.
- The top ranked categories are coloured and labelled; all remaining categories are grey, with rank 1 layered above all other curves.
- Automatic analysis filename suffixes are capped at three parameter tokens.
- Every plot generated through the shared plotting system now receives a metadata TSV sidecar containing the full command invocation and parameters.
- Corrected category-aware aggregate documentation formatting.
