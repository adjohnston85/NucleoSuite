# NucleoSuite 0.8.7

This release adds category-aware regional aggregation.

## Aggregate categories

`nucleosuite aggregate` now accepts `--category-col N`. Each unique value in the selected one-based region-BED column is analysed independently with the existing aggregate pipeline. Each category receives its own aggregate profile, heatmap, processing summary and, by default, independent unified peak calling plus positive- and negative-direction NRL regressions.

A combined category-profile TSV and overlay figure are written after all category analyses complete. A combined category NRL summary reports category, valid region count, directional repeat length, peak count, fit statistics and quality status. Multicontig aggregation remains exact because each category reuses the existing per-contig aggregate/combine implementation.

Combined category output names encode the selected category column together with the central aggregate and NRL parameters, preventing parameter changes from silently overwriting earlier category plots.

For a BED6 TSS file with categories in column 4:

```bash
nucleosuite aggregate \
  --bigwig sample_PNS.bw \
  --region-bed gene_sets_final_tss.bed \
  --category-col 4 \
  --strand-col 6 \
  --missing-strand error \
  --window-half 2500 \
  --output-dir sample_TSS \
  --output-prefix sample_PNS_TSS
```

`--state-bed` help and documentation now explicitly describe it as an inclusion mask rather than a category grouping option.
