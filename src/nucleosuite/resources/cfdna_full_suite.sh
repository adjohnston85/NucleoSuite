#!/usr/bin/env bash
# Run a cfDNA-specific NucleoSuite workflow for paired-end fragment data.

set -Eeuo pipefail
IFS=$'\n\t'

NUCLEOSUITE_BIN="${NUCLEOSUITE_BIN:-nucleosuite}"
PYTHON_BIN="${PYTHON_BIN:-python}"
PLOT_EXT="${NUCLEOSUITE_PLOT_FORMAT:-png}"

BAM_INPUTS=()
BAMS=()
FRAGMENT_INPUTS=()
FRAGMENTS=()
PROVENANCE_BAMS=()
PROVENANCE_FRAGMENTS=()
INPUT_MODE=""
SAMPLE_NAME=""
FASTA=""
BLACKLIST_BED=""
NO_BLACKLIST=0
CTCF_BED=""
STATES_BED=""
GENES_BED=""
GENE_SET_CONFIG=""
EXPRESSION=""
EXPRESSION_VALUE_COLUMN="nTPM"
EXPRESSION_GENE_COLUMN="Gene"
EXPRESSION_NAME_COLUMN="Gene name"
EXPRESSION_PROFILE_COLUMN="Cell line"
EXPRESSION_FOCUS_PROFILES=()
TSS_EXPRESSION_RESOURCE=""
TSS_EXPRESSION_TISSUE="bone_marrow"
TSS_EXPRESSION_WINDOW=2000
RESOURCE_SET=""
VENN_SETS="active_genes,weak_genes,repressed_genes"
OUTDIR=""
CONTIGS=("autosomes")
INTERVAL_FORMAT="both"
ANALYSIS_CHROM_SIZES_SOURCE=""

PNS_FRAG_LOWER=137
PNS_FRAG_UPPER=197
PNS_MODE_LENGTH=167
BIGBED_SCORE_SCALE=1
PNS_SMOOTH_WINDOW=0
PNS_SMOOTH_ORDER=2
PNS_MAX_NEG_RUN=0

EXACT_LENGTHS=(145 161 167)
RANGE_SPECS=("144:146" "160:162" "166:168")
MAX_DUPLICATES=1
MAX_PER_COORDINATE=0
DEDUP_SCOPE="all_bams"
EVEN_DYAD="split"

DAC_DMAX=2000
DAC_WINDOW_SIZE=100000
DAC_ALGORITHM="auto"
NRL_MIN_DISTANCE=1
NRL_MAX_DISTANCE=1500
NRL_PEAK_RESOLUTION=160
SHORT_PERIODICITY_MIN=1
SHORT_PERIODICITY_MAX=144
INTERMEDIATE_PERIODICITY_RESOLUTION=8
INTERMEDIATE_PERIODICITY_MAX=220

AGGREGATE_WINDOW_HALF=2500
STATES_LABEL_COLUMN=4
FRAG_COUNT_MIN=80
FRAG_COUNT_MAX=1000
FRAG_PLOT_MIN=80
FRAG_PLOT_MAX=1000
HEATMAP_MIN_FRAG=80
HEATMAP_MAX_FRAG=500
HEATMAP_NORMALIZATION="fragment-zscore"

CTCF_FLANK=2000
REGION_PEAK_FLANK=2000
DISTANCE_ADJACENT_MAX=500
DISTANCE_LONG_MAX=1500
DISTANCE_LONG_MAX_ORDER=7
STATE_DISTANCE_MAX=500
STATE_DISTANCE_SMOOTH_WINDOW=21
STATE_DISTANCE_SMOOTH_ORDER=2
POSITION_PERCENTILE_INTERVAL=25
SCORE_Z_LIMIT=10
DISTANCE_HISTOGRAM_X_MAX=300
PERCENTILE_BOXPLOT_Y_MAX=500
PEAK_SCORE_NORMALIZATION="count"
POSITIVE_RUNS_THRESHOLD=0
POSITIVE_RUNS_CHUNK_SIZE=1000000
POSITIVE_RUNS_MIN_LENGTH=1
POSITIVE_RUNS_MAX_LENGTH=0
POSITIVE_RUNS_PLOT_X_MAX=550
POSITIVE_RUNS_NORMALIZATION="count"
GENE_FFT_WINDOW=10000
GENE_FFT_PERIOD_MIN=120
GENE_FFT_PERIOD_MAX=280
GENE_FFT_RANKING_PERIODS="193,196,199"
RANDOMIZE_SEED=12345
RANDOMIZE_SEARCH_WINDOW=100000
RANDOMIZE_FALLBACK="uniform"

SKIP_NRL=0
SKIP_FRAGMENT_HEATMAP=0
SKIP_REGION_EXTRACT=0
SKIP_GENE_EXPRESSION=0
SKIP_TSS_EXPRESSION_QUINTILES=0
SKIP_POSITIVE_RUNS=0
SKIP_PEAK_SCORE_FREQUENCY=0
FORCE=0
REUSE_EXISTING_OUTPUTS=0
VALIDATE_ONLY=0
DRY_RUN=0
COMBINE_PREREQUISITES_ONLY=0
WRAPPER_ANALYSIS_SCOPE="combined-only"
RUN_MODE="observed"
RANDOMIZED_INPUT_READY=0
TRUST_EXISTING_OUTPUTS=0

plotting_usage() {
cat <<'EOF'
Shared plot customization:
  --plot-format {png,svg}       Figure format. Default: png.
  --plot-width N                Figure width in inches.
  --plot-height N               Figure height in inches.
  --plot-dpi N                  Figure resolution; relevant mainly to PNG.
  --plot-title TEXT             Override generated plot titles.
  --no-plot-title               Remove plot titles.
  --plot-x-label TEXT           Override x-axis labels.
  --plot-y-label TEXT           Override y-axis labels.
  --plot-font-size N            Base plot font size.
  --plot-grid {none,x,y,both}   Tick-aligned grid lines crossing the plot.
  --plot-grid-color COLOR       Grid-line color.
  --plot-grid-alpha N           Grid-line opacity (0-1).
  --plot-grid-width N           Grid-line width.
  --plot-x-min/--plot-x-max N   Override displayed x-axis limits.
  --plot-y-min/--plot-y-max N   Override displayed y-axis limits.
  --plot-line-width N           Data-line width.
  --plot-line-color COLOR       Single-series line color override.
  --plot-fill-color COLOR       Single-series fill/bar color override.
  --plot-points                 Show point markers on line plots.
  --no-plot-points              Hide point markers.
  --plot-point-size N           Point-marker size.
  --plot-point-fill COLOR       Point-marker fill color.
  --plot-point-edge COLOR       Point-marker outline color.
  --plot-point-edge-width N     Point-marker outline width.
  --plot-point-shape SHAPE      circle, square, triangle or diamond.
  --plot-label-points MODE      none, peaks or all. NRL defaults to peaks; DAC defaults to none.
  --plot-point-label-value V    x, y or both. Default: x.
  --plot-point-label-offset N   Vertical label offset above the point.
  --plot-legend/--no-plot-legend
                                Show or hide legends where available.
  --plot-legend-position POS    best, upper-right, upper-left, lower-right,
                                lower-left or outside-right.
  --plot-x-tick-rotation N      X tick-label rotation in degrees.
  --plot-y-tick-rotation N      Y tick-label rotation in degrees.
  --plot-transparent            Save figures with a transparent background.
EOF
}

usage() {
cat <<'EOF'
Usage:
  nucleosuite cfdna-suite \
      --bam sample.bam \
      --fasta genome.fa \
      --resource-set hg19-gm12878 \
      --outdir results \
      [options]

Inputs and execution:
  --bam FILE_OR_GLOB [MORE ...] Coordinate-sorted paired-end BAM input. Mutually
                                exclusive with --fragments.
  --fragments FILE [MORE ...]   Fragment BED, BED.gz or bigBed input. Mutually
                                exclusive with --bam.
  --sample-name NAME            Output sample name. Default: derived from input filenames.
  --fasta FILE                  Matching indexed reference FASTA.
  --blacklist-bed FILE          Override the assembly-specific blacklist.
  --no-blacklist                Disable blacklist filtering. hg19 v2 is otherwise
                                enabled automatically for hg19/GRCh37 inputs.
  --outdir DIR                  Output directory.
  --cores N                     Process up to N contigs concurrently. Default: 1.
  --combine-cores N             Default worker count for streaming combines unless overridden.
  --streaming-combine-cores N   Memory-light combine workers. Default: --cores.
  --indexed-combine-cores N     BigWig/BigBed combine workers. Default: 1.
  --combine-chunk-bp N          BigWig genomic combine chunk. Default: 100000.
  --analysis-cores N            Memory-light analysis workers. Default: --cores.
  --memory-intensive-analysis-cores N
                                Memory-heavy analysis workers. Default: 1.
                                Neither value may exceed --cores.
  --combine-bigwig-method M      direct streams per-contig BigWigs directly into the
                                combined BigWig; bedgraph writes validated staged
                                bedGraphs during per-contig track generation. Default: direct.

Plot customization:
  --help-plotting               Show the full shared plot customization options and exit.

Multicontig execution:
  --analysis-scope VALUE        combined-only (default) makes per-contig workers generate
                                only combine prerequisites, then runs all analyses once on
                                the combined selected chromosomes. Use
                                per-contig-and-combined to also run downstream analyses for
                                each selected contig.

Resources:
  --resource-set NAME           Bundled resource collection, such as hg19-gm12878.
  --ctcf-bed FILE               CTCF coordinates used for aggregate profiles.
  --states-bed FILE             ChromHMM segmentation used for gene-set generation.
  --genes-bed FILE              Gene regions used for gene-set generation.
  --gene-set-config FILE        Gene-set rule table. The bundled default produces
                                active, weak, repressed and leftover gene sets.
  --expression FILE             Long-format expression TSV. When supplied, expression
                                analysis runs automatically in the combined workflow.
  --expression-value-column N   Expression value column. Default: nTPM.
  --expression-gene-column N    Ensembl gene ID column. Default: Gene.
  --expression-name-column N    Gene-name column. Default: Gene name.
  --expression-profile-column N Expression profile/cell-line column. Default: Cell line.
  --expression-focus-profile N  Profile highlighted in plots; may be repeated.
  --tss-expression-resource F   Tissue-expression TSV/TSV.gz for TSS quintiles. Default: bundled HPA tissue consensus.
  --tss-expression-tissue N     Tissue/profile selector; use underscores for spaces. Default: bone_marrow.
  --tss-expression-window N     Bases on each side of TSS. Default: 2000.
  --venn-sets NAMES             Comma-separated two or three candidate gene sets.
                                Default: active_genes,weak_genes,repressed_genes.
  --states-label-column N       One-based ChromHMM label column. Default: 4.

Fragment classes:
  --score-frag-lower N            PNS lower fragment length. Default: 137.
  --score-frag-upper N            PNS upper fragment length. Default: 197.
  --score-mode-length N           PNS modal fragment length. Default: 167.
  --bigbed-score-scale N        PNS peak multiplier used for integer bigBed scores. Default: 1.
  --score-smooth-window N         Optional Savitzky-Golay window; 0 disables. Default: 0.
  --score-smooth-order N          Optional smoothing polynomial order. Default: 2.
  --score-max-neg-run N           Zero-or-negative bases bridged within PNS peaks. Default: 0.
  --exact-lengths LIST          Comma-separated exact dyad lengths. Default: 145,161,167.
  --range-lengths LIST          Comma-separated lower-upper ranges. Default:
                                144-146,160-162,166-168.
  --max-duplicates N            Identical-fragment copy limit; 0 disables. Default: 1.
  --max-per-coordinate N        Optional dyad/end coordinate cap; 0 disables. Default: 0.
  --dedup-scope VALUE           all_bams or per_bam. Default: all_bams.
  --even-dyad VALUE             split, left or right for dyad tracks. Default: split.

Analysis:
  --contigs VALUE [MORE ...]    Contigs to analyse. Default: autosomes.
  --aggregate-window-half N     Bases on each side of CTCF centres. Default: 2500.
  --dac-dmax N                  Maximum DAC distance. Default: 2000.
  --dac-window-size N           DAC genomic window size. Default: 100000.
  --dac-algorithm VALUE         auto, sparse, or fft. Default: auto.
  --nrl-min-distance N          Long-range NRL regression lower distance. Default: 1.
  --nrl-max-distance N          Long-range NRL regression upper distance. Default: 1500.
  --nrl-peak-resolution N       Long-range NRL peak resolution in bp. Default: 160.
                                Detection smoothing uses resolution/2.5 and local-max smoothing uses resolution/6, snapped down to 10n+1 windows.
  --skip-nrl                    Skip NRL and periodicity analysis of DAC profiles.
  --ctcf-flank N                Region-extraction half-width. Default: 2000.
  --region-peak-flank N         Peak search flank for region extraction. Default: 2000.
  --distance-adjacent-max N     Maximum adjacent-nucleosome distance. Default: 500.
  --distance-long-max N         Maximum distance for 1-7 order NRL analysis. Default: 1500.
  --distance-long-max-order N   Maximum neighbour order for NRL regression. Default: 7.
  --state-distance-max N        Maximum adjacent distance in ChromHMM overlays. Default: 500.
  --state-distance-smooth-window N  State-overlay Savitzky–Golay window. Default: 21.
  --state-distance-smooth-order N   State-overlay polynomial order. Default: 2.
  --peak-score-normalization V  count, fraction, percent, or density. Default: count.
  --positive-runs-threshold N   A base is positive when its score is greater than N. Default: 0.
  --positive-runs-chunk-size N  BigWig scan chunk size. Default: 1000000.
  --positive-runs-min-length N  Minimum retained positive-run length. Default: 1.
  --positive-runs-max-length N  Maximum retained positive-run length; 0 disables. Default: 0.
  --positive-runs-plot-x-max N  Displayed positive-run maximum. Default: 550.
  --positive-runs-normalization V  count, fraction, or percent. Default: count.
  --randomize-seed N            Dinucleotide-matched randomization seed. Default: 12345.
  --randomize-search-window N   Local randomization search window. Default: 100000.
  --randomize-fallback V        uniform or skip. Default: uniform.
  --skip-fragment-heatmap       Skip the fragment-length heatmap.
  --skip-region-extract         Skip CTCF region extraction.
  --skip-tss-expression-quintiles Skip tissue-expression quintile TSS analysis.
  --skip-gene-expression        Skip expression analysis even when --expression is supplied.
  --randomize                   Run a randomized control instead of the observed analysis.
                                The full normal analysis tree is retained and all output
                                names contain _randomized_control.
  --with-randomized-control    Run complete observed and randomized workflows, then
                                append empirical FDR to observed combined peak BEDs.
  --fdr N                       In paired mode, also write BEDs filtered at FDR N.
  --skip-positive-runs          Skip PNS positive-run analysis.
  --skip-peak-score-frequency   Skip peak-score plots for the active run mode.
  --interval-format VALUE       bed, bigbed or both. Default: both.
  --resume                      Reuse only matching completed outputs.
  --force                       Re-run completed steps.
  --dry-run                     Validate inputs and print planned commands without running them.
  -h, --help                    Show this help.

Fragment-length summaries:
  --frag-count-min N            Minimum fragment length counted. Default: 80.
  --frag-count-max N            Maximum fragment length counted. Default: 1000.
  --frag-plot-min N             Fragment-length plot minimum. Default: 80.
  --frag-plot-max N             Plot upper limit; stops at longest count. Default: 1000.
  --heatmap-min-frag N          Heatmap fragment-length minimum. Default: 80.
  --heatmap-max-frag N          Heatmap fragment-length maximum. Default: 500.
  --heatmap-normalization V     fragment-zscore, profile-percent,
                                fragment-percent, profile-minmax,
                                fragment-minmax, or none. Default: fragment-zscore.

Gene-expression FFT settings:
  --gene-fft-window N           Strand-aware window from each TSS. Default: 10000.
  --gene-fft-period-min N       Minimum FFT period. Default: 120.
  --gene-fft-period-max N       Maximum FFT period. Default: 280.
  --gene-fft-ranking-periods N  Comma-separated ranking periods. Default: 193,196,199.
EOF
}

fatal() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "${2:-}" ]] || fatal "$1 requires a value"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) shift; while [[ $# -gt 0 && "$1" != -* ]]; do BAM_INPUTS+=("$1"); shift; done ;;
    --bam=*) BAM_INPUTS+=("${1#*=}"); shift ;;
    --fragments) shift; while [[ $# -gt 0 && "$1" != -* ]]; do FRAGMENT_INPUTS+=("$1"); shift; done ;;
    --fragments=*) FRAGMENT_INPUTS+=("${1#*=}"); shift ;;
    --sample-name) require_value "$@"; SAMPLE_NAME="$2"; shift 2 ;;
    --fasta) require_value "$@"; FASTA="$2"; shift 2 ;;
    --blacklist-bed) require_value "$@"; BLACKLIST_BED="$2"; NO_BLACKLIST=0; shift 2 ;;
    --no-blacklist) NO_BLACKLIST=1; BLACKLIST_BED=""; shift ;;
    --analysis-chrom-sizes-source) require_value "$@"; ANALYSIS_CHROM_SIZES_SOURCE="$2"; shift 2 ;;
    --ctcf-bed) require_value "$@"; CTCF_BED="$2"; shift 2 ;;
    --states-bed) require_value "$@"; STATES_BED="$2"; shift 2 ;;
    --genes-bed) require_value "$@"; GENES_BED="$2"; shift 2 ;;
    --gene-set-config) require_value "$@"; GENE_SET_CONFIG="$2"; shift 2 ;;
    --expression) require_value "$@"; EXPRESSION="$2"; shift 2 ;;
    --expression-value-column) require_value "$@"; EXPRESSION_VALUE_COLUMN="$2"; shift 2 ;;
    --expression-gene-column) require_value "$@"; EXPRESSION_GENE_COLUMN="$2"; shift 2 ;;
    --expression-name-column) require_value "$@"; EXPRESSION_NAME_COLUMN="$2"; shift 2 ;;
    --expression-profile-column) require_value "$@"; EXPRESSION_PROFILE_COLUMN="$2"; shift 2 ;;
    --expression-focus-profile) require_value "$@"; EXPRESSION_FOCUS_PROFILES+=("$2"); shift 2 ;;
    --tss-expression-resource) require_value "$@"; TSS_EXPRESSION_RESOURCE="$2"; shift 2 ;;
    --tss-expression-tissue) require_value "$@"; TSS_EXPRESSION_TISSUE="$2"; shift 2 ;;
    --tss-expression-window) require_value "$@"; TSS_EXPRESSION_WINDOW="$2"; shift 2 ;;
    --venn-sets) require_value "$@"; VENN_SETS="$2"; shift 2 ;;
    --resource-set) require_value "$@"; RESOURCE_SET="$2"; shift 2 ;;
    --outdir) require_value "$@"; OUTDIR="$2"; shift 2 ;;
    --contigs) shift; CONTIGS=(); while [[ $# -gt 0 && "$1" != -* ]]; do CONTIGS+=("$1"); shift; done ;;
    --contigs=*) CONTIGS=("${1#*=}"); shift ;;
    --interval-format) require_value "$@"; INTERVAL_FORMAT="$2"; shift 2 ;;
    --score-frag-lower) require_value "$@"; PNS_FRAG_LOWER="$2"; shift 2 ;;
    --score-frag-upper) require_value "$@"; PNS_FRAG_UPPER="$2"; shift 2 ;;
    --score-mode-length) require_value "$@"; PNS_MODE_LENGTH="$2"; shift 2 ;;
    --bigbed-score-scale) require_value "$@"; BIGBED_SCORE_SCALE="$2"; shift 2 ;;
    --score-smooth-window) require_value "$@"; PNS_SMOOTH_WINDOW="$2"; shift 2 ;;
    --score-smooth-order) require_value "$@"; PNS_SMOOTH_ORDER="$2"; shift 2 ;;
    --score-max-neg-run) require_value "$@"; PNS_MAX_NEG_RUN="$2"; shift 2 ;;
    --exact-lengths) require_value "$@"; IFS=',' read -r -a EXACT_LENGTHS <<< "$2"; shift 2 ;;
    --range-lengths)
      require_value "$@"; IFS=',' read -r -a raw_ranges <<< "$2"; RANGE_SPECS=()
      for value in "${raw_ranges[@]}"; do RANGE_SPECS+=("${value/-/:}"); done; shift 2 ;;
    --max-duplicates) require_value "$@"; MAX_DUPLICATES="$2"; shift 2 ;;
    --max-per-coordinate) require_value "$@"; MAX_PER_COORDINATE="$2"; shift 2 ;;
    --dedup-scope) require_value "$@"; DEDUP_SCOPE="$2"; shift 2 ;;
    --even-dyad) require_value "$@"; EVEN_DYAD="$2"; shift 2 ;;
    --aggregate-window-half) require_value "$@"; AGGREGATE_WINDOW_HALF="$2"; shift 2 ;;
    --states-label-column) require_value "$@"; STATES_LABEL_COLUMN="$2"; shift 2 ;;
    --dac-dmax) require_value "$@"; DAC_DMAX="$2"; shift 2 ;;
    --dac-window-size) require_value "$@"; DAC_WINDOW_SIZE="$2"; shift 2 ;;
    --dac-algorithm) require_value "$@"; DAC_ALGORITHM="$2"; shift 2 ;;
    --nrl-min-distance) require_value "$@"; NRL_MIN_DISTANCE="$2"; shift 2 ;;
    --nrl-max-distance) require_value "$@"; NRL_MAX_DISTANCE="$2"; shift 2 ;;
    --nrl-peak-resolution) require_value "$@"; NRL_PEAK_RESOLUTION="$2"; shift 2 ;;
    --frag-count-min) require_value "$@"; FRAG_COUNT_MIN="$2"; shift 2 ;;
    --frag-count-max) require_value "$@"; FRAG_COUNT_MAX="$2"; shift 2 ;;
    --frag-plot-min) require_value "$@"; FRAG_PLOT_MIN="$2"; shift 2 ;;
    --frag-plot-max) require_value "$@"; FRAG_PLOT_MAX="$2"; shift 2 ;;
    --heatmap-min-frag) require_value "$@"; HEATMAP_MIN_FRAG="$2"; shift 2 ;;
    --heatmap-max-frag) require_value "$@"; HEATMAP_MAX_FRAG="$2"; shift 2 ;;
    --heatmap-normalization) require_value "$@"; HEATMAP_NORMALIZATION="$2"; shift 2 ;;
    --ctcf-flank) require_value "$@"; CTCF_FLANK="$2"; shift 2 ;;
    --region-peak-flank) require_value "$@"; REGION_PEAK_FLANK="$2"; shift 2 ;;
    --distance-adjacent-max) require_value "$@"; DISTANCE_ADJACENT_MAX="$2"; shift 2 ;;
    --distance-long-max) require_value "$@"; DISTANCE_LONG_MAX="$2"; shift 2 ;;
    --distance-long-max-order) require_value "$@"; DISTANCE_LONG_MAX_ORDER="$2"; shift 2 ;;
    --state-distance-max) require_value "$@"; STATE_DISTANCE_MAX="$2"; shift 2 ;;
    --state-distance-smooth-window) require_value "$@"; STATE_DISTANCE_SMOOTH_WINDOW="$2"; shift 2 ;;
    --state-distance-smooth-order) require_value "$@"; STATE_DISTANCE_SMOOTH_ORDER="$2"; shift 2 ;;
    --peak-score-normalization) require_value "$@"; PEAK_SCORE_NORMALIZATION="$2"; shift 2 ;;
    --positive-runs-threshold) require_value "$@"; POSITIVE_RUNS_THRESHOLD="$2"; shift 2 ;;
    --positive-runs-chunk-size) require_value "$@"; POSITIVE_RUNS_CHUNK_SIZE="$2"; shift 2 ;;
    --positive-runs-min-length) require_value "$@"; POSITIVE_RUNS_MIN_LENGTH="$2"; shift 2 ;;
    --positive-runs-max-length) require_value "$@"; POSITIVE_RUNS_MAX_LENGTH="$2"; shift 2 ;;
    --positive-runs-plot-x-max) require_value "$@"; POSITIVE_RUNS_PLOT_X_MAX="$2"; shift 2 ;;
    --positive-runs-normalization) require_value "$@"; POSITIVE_RUNS_NORMALIZATION="$2"; shift 2 ;;
    --gene-fft-window) require_value "$@"; GENE_FFT_WINDOW="$2"; shift 2 ;;
    --gene-fft-period-min) require_value "$@"; GENE_FFT_PERIOD_MIN="$2"; shift 2 ;;
    --gene-fft-period-max) require_value "$@"; GENE_FFT_PERIOD_MAX="$2"; shift 2 ;;
    --gene-fft-ranking-periods) require_value "$@"; GENE_FFT_RANKING_PERIODS="$2"; shift 2 ;;
    --randomize-seed) require_value "$@"; RANDOMIZE_SEED="$2"; shift 2 ;;
    --randomize-search-window) require_value "$@"; RANDOMIZE_SEARCH_WINDOW="$2"; shift 2 ;;
    --randomize-fallback) require_value "$@"; RANDOMIZE_FALLBACK="$2"; shift 2 ;;
    --skip-nrl) SKIP_NRL=1; shift ;;
    --skip-fragment-heatmap) SKIP_FRAGMENT_HEATMAP=1; shift ;;
    --skip-region-extract) SKIP_REGION_EXTRACT=1; shift ;;
    --skip-gene-expression) SKIP_GENE_EXPRESSION=1; shift ;;
    --skip-tss-expression-quintiles) SKIP_TSS_EXPRESSION_QUINTILES=1; shift ;;
    --randomize) RUN_MODE="randomized"; shift ;;
    --randomized-control-input) RUN_MODE="randomized"; RANDOMIZED_INPUT_READY=1; shift ;;
    --provenance-bam) require_value "$@"; PROVENANCE_BAMS+=("$2"); shift 2 ;;
    --provenance-fragment) require_value "$@"; PROVENANCE_FRAGMENTS+=("$2"); shift 2 ;;
    --trusted-combined-prerequisites) TRUST_EXISTING_OUTPUTS=1; shift ;;
    --skip-positive-runs) SKIP_POSITIVE_RUNS=1; shift ;;
    --skip-peak-score-frequency) SKIP_PEAK_SCORE_FREQUENCY=1; shift ;;
    --resume) REUSE_EXISTING_OUTPUTS=1; shift ;;
    --force) FORCE=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    --validate-only) VALIDATE_ONLY=1; shift ;;
    --combine-bigwig-method)
      require_value "$@"
      [[ "$2" == "direct" || "$2" == "bedgraph" || "$2" == "bedgraphs" ]] || fatal "--combine-bigwig-method must be direct or bedgraph"
      shift 2
      ;;
    --combine-bigwig-method=*)
      WRAPPER_BIGWIG_METHOD="${1#*=}"
      [[ "$WRAPPER_BIGWIG_METHOD" == "direct" || "$WRAPPER_BIGWIG_METHOD" == "bedgraph" || "$WRAPPER_BIGWIG_METHOD" == "bedgraphs" ]] || fatal "--combine-bigwig-method must be direct or bedgraph"
      shift
      ;;
    --analysis-scope)
      require_value "$@"
      WRAPPER_ANALYSIS_SCOPE="$2"
      [[ "$WRAPPER_ANALYSIS_SCOPE" == "combined-only" || "$WRAPPER_ANALYSIS_SCOPE" == "per-contig-and-combined" ]] || fatal "invalid --analysis-scope"
      shift 2
      ;;
    --analysis-scope=*)
      WRAPPER_ANALYSIS_SCOPE="${1#*=}"
      [[ "$WRAPPER_ANALYSIS_SCOPE" == "combined-only" || "$WRAPPER_ANALYSIS_SCOPE" == "per-contig-and-combined" ]] || fatal "invalid --analysis-scope"
      shift
      ;;
    --combine-prerequisites-only) COMBINE_PREREQUISITES_ONLY=1; shift ;;
    --cores) shift 2 ;; # consumed by the Python wrapper
    --cores=*) shift ;;
    -h|--help) usage; exit 0 ;;
    --help-plotting) usage; plotting_usage; exit 0 ;;
    *) fatal "Unknown option: $1" ;;
  esac
done

[[ -n "$OUTDIR" ]] || fatal "--outdir is required"
[[ -n "$FASTA" ]] || fatal "--fasta is required"
[[ -f "$FASTA" ]] || fatal "FASTA not found: $FASTA"
[[ -z "$BLACKLIST_BED" || -f "$BLACKLIST_BED" ]] || fatal "Blacklist BED not found: $BLACKLIST_BED"
if [[ ${#BAM_INPUTS[@]} -gt 0 && ${#FRAGMENT_INPUTS[@]} -gt 0 ]]; then fatal "Use either --bam or --fragments, not both"; fi
if [[ ${#BAM_INPUTS[@]} -eq 0 && ${#FRAGMENT_INPUTS[@]} -eq 0 ]]; then fatal "Provide --bam or --fragments"; fi
[[ "$MAX_DUPLICATES" =~ ^[0-9]+$ ]] || fatal "--max-duplicates must be a non-negative integer"
[[ "$MAX_PER_COORDINATE" =~ ^[0-9]+$ ]] || fatal "--max-per-coordinate must be a non-negative integer"
[[ "$DEDUP_SCOPE" == "all_bams" || "$DEDUP_SCOPE" == "per_bam" ]] || fatal "invalid --dedup-scope"
[[ "$EVEN_DYAD" == "split" || "$EVEN_DYAD" == "left" || "$EVEN_DYAD" == "right" ]] || fatal "invalid --even-dyad"
[[ "$INTERVAL_FORMAT" =~ ^(bed|bigbed|both)$ ]] || fatal "invalid --interval-format"
[[ "$PNS_SMOOTH_WINDOW" =~ ^[0-9]+$ ]] || fatal "--score-smooth-window must be a non-negative integer"
[[ "$PNS_SMOOTH_ORDER" =~ ^[0-9]+$ ]] || fatal "--score-smooth-order must be a non-negative integer"
[[ "$PNS_MAX_NEG_RUN" =~ ^[0-9]+$ ]] || fatal "--score-max-neg-run must be a non-negative integer"
if [[ "$PNS_SMOOTH_WINDOW" -gt 0 ]]; then
  [[ "$PNS_SMOOTH_WINDOW" -ge 3 && $((PNS_SMOOTH_WINDOW % 2)) -eq 1 ]] || fatal "--score-smooth-window must be 0 or an odd integer of at least 3"
  [[ "$PNS_SMOOTH_ORDER" -lt "$PNS_SMOOTH_WINDOW" ]] || fatal "--score-smooth-order must be smaller than --score-smooth-window"
fi
for pair in \
  "--ctcf-flank:$CTCF_FLANK" "--region-peak-flank:$REGION_PEAK_FLANK" \
  "--distance-adjacent-max:$DISTANCE_ADJACENT_MAX" "--distance-long-max:$DISTANCE_LONG_MAX" "--distance-long-max-order:$DISTANCE_LONG_MAX_ORDER" \
  "--state-distance-max:$STATE_DISTANCE_MAX" \
  "--positive-runs-chunk-size:$POSITIVE_RUNS_CHUNK_SIZE" \
  "--positive-runs-min-length:$POSITIVE_RUNS_MIN_LENGTH" "--positive-runs-max-length:$POSITIVE_RUNS_MAX_LENGTH" \
  "--positive-runs-plot-x-max:$POSITIVE_RUNS_PLOT_X_MAX" "--gene-fft-window:$GENE_FFT_WINDOW" \
  "--gene-fft-period-min:$GENE_FFT_PERIOD_MIN" "--gene-fft-period-max:$GENE_FFT_PERIOD_MAX" \
  "--randomize-seed:$RANDOMIZE_SEED" "--randomize-search-window:$RANDOMIZE_SEARCH_WINDOW"; do
  [[ "${pair#*:}" =~ ^[0-9]+$ ]] || fatal "${pair%%:*} must be a non-negative integer"
done
[[ "$PEAK_SCORE_NORMALIZATION" =~ ^(count|fraction|percent|density)$ ]] || fatal "invalid --peak-score-normalization"
[[ "$POSITIVE_RUNS_NORMALIZATION" =~ ^(count|fraction|percent)$ ]] || fatal "invalid --positive-runs-normalization"
[[ "$RANDOMIZE_FALLBACK" =~ ^(uniform|skip)$ ]] || fatal "invalid --randomize-fallback"
[[ "$STATE_DISTANCE_SMOOTH_WINDOW" =~ ^[0-9]+$ && $((STATE_DISTANCE_SMOOTH_WINDOW % 2)) -eq 1 ]] || fatal "--state-distance-smooth-window must be odd"
[[ "$STATE_DISTANCE_SMOOTH_ORDER" =~ ^[0-9]+$ && "$STATE_DISTANCE_SMOOTH_ORDER" -lt "$STATE_DISTANCE_SMOOTH_WINDOW" ]] || fatal "invalid --state-distance-smooth-order"
[[ "$GENE_FFT_PERIOD_MIN" -lt "$GENE_FFT_PERIOD_MAX" ]] || fatal "gene FFT minimum period must be less than maximum"

expand_inputs() {
  local token match
  if [[ ${#BAM_INPUTS[@]} -gt 0 ]]; then
    for token in "${BAM_INPUTS[@]}"; do
      if [[ -f "$token" ]]; then BAMS+=("$(realpath "$token")");
      else while IFS= read -r match; do BAMS+=("$(realpath "$match")"); done < <(compgen -G "$token" || true); fi
    done
    [[ ${#BAMS[@]} -gt 0 ]] || fatal "No BAM files matched"
    mapfile -t BAMS < <(printf '%s\n' "${BAMS[@]}" | sort -V -u)
    INPUT_MODE="bam"
  else
    for token in "${FRAGMENT_INPUTS[@]}"; do
      if [[ -f "$token" ]]; then FRAGMENTS+=("$(realpath "$token")");
      else while IFS= read -r match; do FRAGMENTS+=("$(realpath "$match")"); done < <(compgen -G "$token" || true); fi
    done
    [[ ${#FRAGMENTS[@]} -gt 0 ]] || fatal "No fragment files matched"
    mapfile -t FRAGMENTS < <(printf '%s\n' "${FRAGMENTS[@]}" | sort -V -u)
    INPUT_MODE="fragments"
  fi
}
expand_inputs
if [[ "${#PROVENANCE_BAMS[@]}" -eq 0 && "${#PROVENANCE_FRAGMENTS[@]}" -eq 0 ]]; then
  if [[ "$INPUT_MODE" == "bam" ]]; then
    PROVENANCE_BAMS=("${BAMS[@]}")
  else
    PROVENANCE_FRAGMENTS=("${FRAGMENTS[@]}")
  fi
fi

if [[ -n "$RESOURCE_SET" ]]; then
  "$NUCLEOSUITE_BIN" resources validate --resource-set "$RESOURCE_SET" >/dev/null
  [[ -n "$CTCF_BED" ]] || CTCF_BED="$($NUCLEOSUITE_BIN resources show ctcf --resource-set "$RESOURCE_SET")"
  [[ -n "$STATES_BED" ]] || STATES_BED="$($NUCLEOSUITE_BIN resources show states --resource-set "$RESOURCE_SET")"
  [[ -n "$GENES_BED" ]] || GENES_BED="$($NUCLEOSUITE_BIN resources show genes --resource-set "$RESOURCE_SET")"
  [[ -n "$GENE_SET_CONFIG" ]] || GENE_SET_CONFIG="$($NUCLEOSUITE_BIN resources show gene_set_config --resource-set "$RESOURCE_SET")"
  [[ -n "$TSS_EXPRESSION_RESOURCE" ]] || TSS_EXPRESSION_RESOURCE="$($NUCLEOSUITE_BIN resources show tissue_expression --resource-set "$RESOURCE_SET")"
fi
[[ -n "$TSS_EXPRESSION_RESOURCE" ]] || TSS_EXPRESSION_RESOURCE="$($NUCLEOSUITE_BIN resources path hpa-tissue-expression)"
[[ -n "$CTCF_BED" ]] || fatal "--ctcf-bed is required unless --resource-set supplies it"
[[ -f "$CTCF_BED" ]] || fatal "CTCF BED not found: $CTCF_BED"
[[ -z "$STATES_BED" || -f "$STATES_BED" ]] || fatal "states BED not found: $STATES_BED"
[[ -z "$GENES_BED" || -f "$GENES_BED" ]] || fatal "genes BED not found: $GENES_BED"
[[ -z "$GENE_SET_CONFIG" || -f "$GENE_SET_CONFIG" ]] || fatal "gene-set config not found: $GENE_SET_CONFIG"
[[ -z "$EXPRESSION" || -f "$EXPRESSION" ]] || fatal "expression TSV not found: $EXPRESSION"
[[ -f "$TSS_EXPRESSION_RESOURCE" ]] || fatal "TSS expression resource not found: $TSS_EXPRESSION_RESOURCE"
[[ "$TSS_EXPRESSION_WINDOW" =~ ^[1-9][0-9]*$ ]] || fatal "--tss-expression-window must be a positive integer"
if [[ -n "$EXPRESSION" && -z "$GENES_BED" ]]; then fatal "--expression requires --genes-bed or a resource set that supplies genes"; fi
if [[ -n "$STATES_BED" && -n "$GENES_BED" && -z "$GENE_SET_CONFIG" ]]; then
  GENE_SET_CONFIG="$($NUCLEOSUITE_BIN resources path default-gene-sets)"
fi

if [[ -z "$SAMPLE_NAME" ]]; then
  if [[ "$INPUT_MODE" == "bam" ]]; then
    if [[ ${#BAMS[@]} -eq 1 ]]; then SAMPLE_NAME="$(basename "${BAMS[0]}" .bam)"; else SAMPLE_NAME="multi_bam"; fi
  else
    if [[ ${#FRAGMENTS[@]} -eq 1 ]]; then SAMPLE_NAME="$(basename "${FRAGMENTS[0]}")"; else SAMPLE_NAME="multi_fragments"; fi
  fi
fi
SAMPLE_NAME="${SAMPLE_NAME//[^A-Za-z0-9._-]/_}"
SAMPLE="$SAMPLE_NAME"
if [[ "$RUN_MODE" == "randomized" && "$SAMPLE" != *_randomized_control* ]]; then
  SAMPLE="${SAMPLE}_randomized_control"
  SAMPLE_NAME="$SAMPLE"
fi
SUPPORT_PREFIX=""
if [[ "$RUN_MODE" == "randomized" ]]; then
  SUPPORT_PREFIX="${SAMPLE}_"
fi

BLACKLIST_ARGS=()
[[ -n "$BLACKLIST_BED" ]] && BLACKLIST_ARGS=(--blacklist-bed "$BLACKLIST_BED")
ANALYSIS_INPUT_ARGS=()
SOURCE_REFERENCE_ARGS=()
if [[ "$INPUT_MODE" == "bam" ]]; then ANALYSIS_INPUT_ARGS=(-b "${BAMS[@]}");
else ANALYSIS_INPUT_ARGS=(--fragments "${FRAGMENTS[@]}"); SOURCE_REFERENCE_ARGS=(--fasta "$FASTA"); fi
SOURCE_INPUT_MODE="$INPUT_MODE"
SOURCE_INPUT_ARGS=("${ANALYSIS_INPUT_ARGS[@]}")

if [[ "$VALIDATE_ONLY" -eq 1 || "$DRY_RUN" -eq 1 ]]; then
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf 'NucleoSuite cfDNA suite dry run\nmode\t%s\noutdir\t%s\nblacklist\t%s\n' \
      "$RUN_MODE" "$OUTDIR" "${BLACKLIST_BED:-auto-if-hg19}"
    printf 'stages\tsetup,tracks,scaling,dac,nrl,aggregates,distances,region-extract,fragment-lengths,heatmaps,gene-expression,positive-runs,peak-analysis\n'
  fi
  exit 0
fi

mkdir -p "$OUTDIR"
SETUP_DIR="$OUTDIR/00_setup"
GENE_SET_DIR="$OUTDIR/00_gene_sets"
COMBINED_TRACK_DIR="$OUTDIR/01_combined_tracks"
PNS_DIR="$COMBINED_TRACK_DIR/pns"
SCALED_DIR="$COMBINED_TRACK_DIR/scaled"
EXACT_DIR="$COMBINED_TRACK_DIR/dyads/exact"
RANGE_DIR="$COMBINED_TRACK_DIR/dyads/ranges"
END_EXACT_DIR="$COMBINED_TRACK_DIR/fragment_ends/exact"
END_RANGE_DIR="$COMBINED_TRACK_DIR/fragment_ends/ranges"
SEQ_DIR="$COMBINED_TRACK_DIR/sequence"
DAC_DIR="$OUTDIR/02_dac"
NRL_DIR="$OUTDIR/04_nrl"
AGG_DIR="$OUTDIR/05_ctcf_aggregation"
TSS_AGG_DIR="$OUTDIR/06_tss_aggregation"
TSS_EXPRESSION_DIR="$OUTDIR/06_tss_expression_quintiles"
DIST_DIR="$OUTDIR/07_distances"
REGION_DIR="$OUTDIR/08_region_extract"
FRAG_DIR="$OUTDIR/09_fragment_lengths"
HEATMAP_DIR="$OUTDIR/10_fragment_heatmaps"
GENE_EXPRESSION_DIR="$OUTDIR/11_gene_expression"
POSITIVE_RUNS_DIR="$OUTDIR/12_positive_runs"
PEAK_ANALYSIS_DIR="$OUTDIR/13_peak_analysis"
PEAK_SCORE_DIR="$PEAK_ANALYSIS_DIR/score_frequencies"
LOG_DIR="$OUTDIR/logs"
mkdir -p "$SETUP_DIR" "$GENE_SET_DIR" "$COMBINED_TRACK_DIR" "$PNS_DIR" "$SCALED_DIR"   "$EXACT_DIR" "$RANGE_DIR" "$END_EXACT_DIR" "$END_RANGE_DIR" "$SEQ_DIR" "$DAC_DIR"   "$NRL_DIR" "$AGG_DIR" "$TSS_AGG_DIR" "$DIST_DIR" "$REGION_DIR" "$FRAG_DIR" "$HEATMAP_DIR"   "$GENE_EXPRESSION_DIR" "$POSITIVE_RUNS_DIR" "$PEAK_ANALYSIS_DIR" "$PEAK_SCORE_DIR" "$LOG_DIR"
DONE_DIR="$OUTDIR/.done"
mkdir -p "$DONE_DIR"

PARAMETERS="$SETUP_DIR/${SUPPORT_PREFIX}run_parameters.tsv"
{
  printf 'parameter\tvalue\n'
  for name in SAMPLE_NAME RUN_MODE FASTA BLACKLIST_BED NO_BLACKLIST CTCF_BED STATES_BED GENES_BED GENE_SET_CONFIG \
    EXPRESSION EXPRESSION_VALUE_COLUMN EXPRESSION_GENE_COLUMN EXPRESSION_NAME_COLUMN EXPRESSION_PROFILE_COLUMN \
    TSS_EXPRESSION_RESOURCE TSS_EXPRESSION_TISSUE TSS_EXPRESSION_WINDOW RESOURCE_SET \
    VENN_SETS OUTDIR INTERVAL_FORMAT PNS_FRAG_LOWER PNS_FRAG_UPPER PNS_MODE_LENGTH PNS_SMOOTH_WINDOW PNS_SMOOTH_ORDER PNS_MAX_NEG_RUN \
    MAX_DUPLICATES MAX_PER_COORDINATE DEDUP_SCOPE EVEN_DYAD DAC_DMAX DAC_WINDOW_SIZE DAC_ALGORITHM NRL_MIN_DISTANCE NRL_MAX_DISTANCE NRL_PEAK_RESOLUTION \
    AGGREGATE_WINDOW_HALF STATES_LABEL_COLUMN FRAG_COUNT_MIN FRAG_COUNT_MAX FRAG_PLOT_MIN \
    FRAG_PLOT_MAX HEATMAP_MIN_FRAG HEATMAP_MAX_FRAG HEATMAP_NORMALIZATION CTCF_FLANK REGION_PEAK_FLANK \
    DISTANCE_ADJACENT_MAX DISTANCE_LONG_MAX DISTANCE_LONG_MAX_ORDER STATE_DISTANCE_MAX STATE_DISTANCE_SMOOTH_WINDOW STATE_DISTANCE_SMOOTH_ORDER \
    PEAK_SCORE_NORMALIZATION POSITIVE_RUNS_THRESHOLD POSITIVE_RUNS_CHUNK_SIZE POSITIVE_RUNS_MIN_LENGTH \
    POSITIVE_RUNS_MAX_LENGTH POSITIVE_RUNS_PLOT_X_MAX POSITIVE_RUNS_NORMALIZATION GENE_FFT_WINDOW GENE_FFT_PERIOD_MIN \
    GENE_FFT_PERIOD_MAX GENE_FFT_RANKING_PERIODS RANDOMIZE_SEED RANDOMIZE_SEARCH_WINDOW RANDOMIZE_FALLBACK SKIP_NRL SKIP_FRAGMENT_HEATMAP SKIP_REGION_EXTRACT SKIP_GENE_EXPRESSION \
    SKIP_TSS_EXPRESSION_QUINTILES SKIP_POSITIVE_RUNS SKIP_PEAK_SCORE_FREQUENCY COMBINE_PREREQUISITES_ONLY; do
    printf '%s\t%s\n' "$name" "${!name}"
  done
  for index in "${!BAMS[@]}"; do printf 'BAM_%s\t%s\n' "$((index + 1))" "${BAMS[$index]}"; done
  for index in "${!FRAGMENTS[@]}"; do printf 'FRAGMENT_%s\t%s\n' "$((index + 1))" "${FRAGMENTS[$index]}"; done
  printf 'SOURCE_BAM_COUNT\t%s\n' "${#PROVENANCE_BAMS[@]}"
  for index in "${!PROVENANCE_BAMS[@]}"; do printf 'SOURCE_BAM_%s\t%s\n' "$((index + 1))" "${PROVENANCE_BAMS[$index]}"; done
  printf 'SOURCE_FRAGMENT_COUNT\t%s\n' "${#PROVENANCE_FRAGMENTS[@]}"
  for index in "${!PROVENANCE_FRAGMENTS[@]}"; do printf 'SOURCE_FRAGMENT_%s\t%s\n' "$((index + 1))" "${PROVENANCE_FRAGMENTS[$index]}"; done
  for index in "${!EXACT_LENGTHS[@]}"; do printf 'EXACT_LENGTH_%s\t%s\n' "$((index + 1))" "${EXACT_LENGTHS[$index]}"; done
  for index in "${!RANGE_SPECS[@]}"; do printf 'RANGE_SPEC_%s\t%s\n' "$((index + 1))" "${RANGE_SPECS[$index]}"; done
  for index in "${!EXPRESSION_FOCUS_PROFILES[@]}"; do printf 'EXPRESSION_FOCUS_PROFILE_%s\t%s\n' "$((index + 1))" "${EXPRESSION_FOCUS_PROFILES[$index]}"; done
  printf 'CONTIGS\t%s\n' "$(IFS=,; echo "${CONTIGS[*]}")"
} > "$PARAMETERS"

append_input_identities() {
  local parameter_file="$1"; shift
  [[ "$#" -gt 0 ]] || return 0
  "$PYTHON_BIN" - "$parameter_file" "$@" <<'PYIDENTITY'
import hashlib, json, sys
from pathlib import Path

parameter_path = Path(sys.argv[1])
items = sys.argv[2:]
if len(items) % 2:
    raise SystemExit("Input identity arguments must be LABEL PATH pairs")

def fingerprint(path_text):
    path = Path(path_text)
    if not path.is_file():
        raise SystemExit(f"Provenance input is not a readable file: {path}")
    stat = path.stat()
    size = stat.st_size
    digest = hashlib.sha256()
    block = 1024 * 1024
    if size <= 16 * block:
        method = "sha256_full"
        with path.open("rb") as handle:
            while True:
                chunk = handle.read(block)
                if not chunk:
                    break
                digest.update(chunk)
    else:
        method = "sha256_sampled_first_middle_last_1MiB"
        with path.open("rb") as handle:
            for offset in (0, max(0, size // 2 - block // 2), max(0, size - block)):
                handle.seek(offset)
                digest.update(offset.to_bytes(8, "little"))
                digest.update(handle.read(block))
        digest.update(size.to_bytes(16, "little"))
    return {
        "path": str(path.resolve()),
        "size_bytes": size,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
        "device": stat.st_dev,
        "inode": stat.st_ino,
        "fingerprint_method": method,
        "sha256": digest.hexdigest(),
    }

with parameter_path.open("a", encoding="utf-8") as handle:
    for index in range(0, len(items), 2):
        label, path = items[index], items[index + 1]
        if path:
            handle.write(
                f"INPUT_IDENTITY_{label}\t"
                + json.dumps(fingerprint(path), sort_keys=True, separators=(",", ":"))
                + "\n"
            )
PYIDENTITY
}

INPUT_IDENTITY_ARGS=(FASTA "$FASTA" CTCF "$CTCF_BED" TSS_EXPRESSION "$TSS_EXPRESSION_RESOURCE")
[[ -z "$STATES_BED" ]] || INPUT_IDENTITY_ARGS+=(STATES "$STATES_BED")
[[ -z "$GENES_BED" ]] || INPUT_IDENTITY_ARGS+=(GENES "$GENES_BED")
[[ -z "$GENE_SET_CONFIG" ]] || INPUT_IDENTITY_ARGS+=(GENE_SET_CONFIG "$GENE_SET_CONFIG")
[[ -z "$EXPRESSION" ]] || INPUT_IDENTITY_ARGS+=(EXPRESSION "$EXPRESSION")
[[ -z "$ANALYSIS_CHROM_SIZES_SOURCE" ]] || INPUT_IDENTITY_ARGS+=(ANALYSIS_CHROM_SIZES_SOURCE "$ANALYSIS_CHROM_SIZES_SOURCE")
for index in "${!BAMS[@]}"; do INPUT_IDENTITY_ARGS+=("BAM_$((index + 1))" "${BAMS[$index]}"); done
for index in "${!FRAGMENTS[@]}"; do INPUT_IDENTITY_ARGS+=("FRAGMENT_$((index + 1))" "${FRAGMENTS[$index]}"); done
for index in "${!PROVENANCE_BAMS[@]}"; do INPUT_IDENTITY_ARGS+=("SOURCE_BAM_$((index + 1))" "${PROVENANCE_BAMS[$index]}"); done
for index in "${!PROVENANCE_FRAGMENTS[@]}"; do INPUT_IDENTITY_ARGS+=("SOURCE_FRAGMENT_$((index + 1))" "${PROVENANCE_FRAGMENTS[$index]}"); done
append_input_identities "$PARAMETERS" "${INPUT_IDENTITY_ARGS[@]}"
PARAM_HASH="$($PYTHON_BIN - "$PARAMETERS" <<'PY'
import hashlib, sys
print(hashlib.sha256(open(sys.argv[1], 'rb').read()).hexdigest()[:16])
PY
)"

INTERVAL_EXT="bed"
[[ "$INTERVAL_FORMAT" == "bigbed" ]] && INTERVAL_EXT="bb"

PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
FAILED_STEPS=()
FAILED_STEP_STATUS=()
FAILED_STEP_REASON=()
FAILED_STEP_LOG=()

expected_exists() {
  local pattern="$1"
  compgen -G "$pattern" >/dev/null 2>&1
}
record_failed_step() {
  local name="$1" status="$2" reason="$3" log="$4"
  FAIL_COUNT=$((FAIL_COUNT + 1))
  FAILED_STEPS+=("$name")
  FAILED_STEP_STATUS+=("$status")
  FAILED_STEP_REASON+=("$reason")
  FAILED_STEP_LOG+=("$log")
}
run_step() {
  local name="$1" expected="$2"; shift 2
  local log="$LOG_DIR/${SUPPORT_PREFIX}${name}.log" status=0 marker="$DONE_DIR/${SUPPORT_PREFIX}${PARAM_HASH}_${name}.done"
  if [[ "$FORCE" -eq 0 && "$REUSE_EXISTING_OUTPUTS" -eq 1 && ( -f "$marker" || "$TRUST_EXISTING_OUTPUTS" -eq 1 ) ]] && \
     { [[ -e "$expected" ]] || expected_exists "$expected"; }; then
    echo "[SKIP] $name"
    SKIP_COUNT=$((SKIP_COUNT + 1))
    return 0
  fi
  rm -f "$marker"
  echo "[RUN ] $name"
  { printf '[CMD ]'; printf ' %q' "$@"; printf '\n'; } >"$log"
  if "$@" >>"$log" 2>&1; then
    status=0
  else
    status=$?
  fi
  if [[ "$status" -ne 0 ]]; then
    echo "[FAIL] $name exited with status $status (see $log)" >&2
    record_failed_step "$name" "$status" "command_failed" "$log"
    return 0
  fi
  if [[ ! -e "$expected" ]] && ! expected_exists "$expected"; then
    echo "[FAIL] $name did not create $expected (see $log)" >&2
    record_failed_step "$name" "0" "missing_expected_output" "$log"
    return 0
  fi
  PASS_COUNT=$((PASS_COUNT + 1))
  touch "$marker"
  echo "[PASS] $name"
  return 0
}

# Combined-only analytical stages can contain hundreds of independent jobs.
# The Python wrapper supplies one global worker budget through this variable.
ASYNC_CORES="${NUCLEOSUITE_SUITE_CORES:-1}"
[[ "$ASYNC_CORES" =~ ^[1-9][0-9]*$ ]] || ASYNC_CORES=1
MEMORY_ASYNC_CORES="${NUCLEOSUITE_MEMORY_ANALYSIS_CORES:-1}"
[[ "$MEMORY_ASYNC_CORES" =~ ^[1-9][0-9]*$ ]] || MEMORY_ASYNC_CORES=1
ASYNC_STATUS_DIR="$LOG_DIR/parallel_status"
mkdir -p "$ASYNC_STATUS_DIR"
ASYNC_NAMES=()
ASYNC_EXPECTED=()
ASYNC_LOGS=()
ASYNC_STATUS_FILES=()
ASYNC_PIDS=()

finish_queued_step() {
  local index="$1" name log status_file status reason marker
  name="${ASYNC_NAMES[$index]}"
  log="${ASYNC_LOGS[$index]}"
  status_file="${ASYNC_STATUS_FILES[$index]}"
  marker="$DONE_DIR/${SUPPORT_PREFIX}${PARAM_HASH}_${name}.done"
  status=2; reason="missing_parallel_status"
  if [[ -s "$status_file" ]]; then
    IFS=$'\t' read -r status reason < "$status_file" || true
  fi
  if [[ "$status" == "0" ]]; then
    touch "$marker"
    PASS_COUNT=$((PASS_COUNT + 1))
    echo "[PASS] $name"
  else
    echo "[FAIL] $name exited with status $status (see $log)" >&2
    record_failed_step "$name" "$status" "$reason" "$log"
  fi
  rm -f "$status_file"
  unset 'ASYNC_NAMES[index]' 'ASYNC_EXPECTED[index]' 'ASYNC_LOGS[index]' \
    'ASYNC_STATUS_FILES[index]' 'ASYNC_PIDS[index]'
}

wait_one_queued_step() {
  [[ "${#ASYNC_PIDS[@]}" -gt 0 ]] || return 0
  local completed_pid="" index
  wait -n -p completed_pid "${ASYNC_PIDS[@]}" || true
  for index in "${!ASYNC_PIDS[@]}"; do
    if [[ "${ASYNC_PIDS[$index]}" == "$completed_pid" ]]; then
      finish_queued_step "$index"
      return 0
    fi
  done
  index="${!ASYNC_PIDS[@]}"; index="${index%% *}"
  wait "${ASYNC_PIDS[$index]}" || true
  finish_queued_step "$index"
}

wait_queued_steps() {
  while [[ "${#ASYNC_PIDS[@]}" -gt 0 ]]; do
    wait_one_queued_step
  done
}

queue_step() {
  local name="$1" expected="$2"; shift 2
  if [[ "$ASYNC_CORES" -le 1 ]]; then
    run_step "$name" "$expected" "$@"
    return 0
  fi
  local log="$LOG_DIR/${SUPPORT_PREFIX}${name}.log" marker="$DONE_DIR/${SUPPORT_PREFIX}${PARAM_HASH}_${name}.done"
  local safe_name="${name//[^A-Za-z0-9._-]/_}"
  local status_file="$ASYNC_STATUS_DIR/${SUPPORT_PREFIX}${safe_name}.$$.status"
  if [[ "$FORCE" -eq 0 && "$REUSE_EXISTING_OUTPUTS" -eq 1 && ( -f "$marker" || "$TRUST_EXISTING_OUTPUTS" -eq 1 ) ]] && \
     { [[ -e "$expected" ]] || expected_exists "$expected"; }; then
    echo "[SKIP] $name"
    SKIP_COUNT=$((SKIP_COUNT + 1))
    return 0
  else
    rm -f "$marker"
    echo "[RUN ] $name"
    { printf '[CMD ]'; printf ' %q' "$@"; printf '\n'; } >"$log"
    (
      local worker_status=0 worker_reason="command_failed"
      if "$@" >>"$log" 2>&1; then
        worker_status=0
      else
        worker_status=$?
      fi
      if [[ "$worker_status" -eq 0 ]]; then
        if [[ -e "$expected" ]] || expected_exists "$expected"; then
          worker_reason="pass"
        else
          worker_status=2
          worker_reason="missing_expected_output"
        fi
      fi
      printf '%s\t%s\n' "$worker_status" "$worker_reason" > "$status_file"
    ) &
  fi
  ASYNC_NAMES+=("$name")
  ASYNC_EXPECTED+=("$expected")
  ASYNC_LOGS+=("$log")
  ASYNC_STATUS_FILES+=("$status_file")
  ASYNC_PIDS+=("$!")
  if [[ "${#ASYNC_PIDS[@]}" -ge "$ASYNC_CORES" ]]; then
    wait_one_queued_step
  fi
}

queue_memory_step() {
  local normal_cores="$ASYNC_CORES"
  ASYNC_CORES="$MEMORY_ASYNC_CORES"
  queue_step "$@"
  ASYNC_CORES="$normal_cores"
}

# Build analysis chromosome sizes during setup, independently of track or
# regional-resource generation. All downstream stages use this validated file.
ALL_CHROM_SIZES="$SETUP_DIR/${SUPPORT_PREFIX}analysis.chrom.sizes"
CHROM_SIZES="$SETUP_DIR/${SUPPORT_PREFIX}selected.chrom.sizes"
CHROM_SIZE_INPUTS=()
if [[ -n "$ANALYSIS_CHROM_SIZES_SOURCE" ]]; then
    CHROM_SIZE_INPUTS=("$ANALYSIS_CHROM_SIZES_SOURCE")
elif [[ "$INPUT_MODE" == "bam" ]]; then
    CHROM_SIZE_INPUTS=("${BAMS[@]}")
else
    CHROM_SIZE_INPUTS=("${FRAGMENTS[@]}")
fi
run_step "00_chrom_sizes" "$CHROM_SIZES" "$PYTHON_BIN" - \
    "$ANALYSIS_CHROM_SIZES_SOURCE" "$FASTA" "$ALL_CHROM_SIZES" "$CHROM_SIZES" \
    "${CONTIGS[*]}" "$INPUT_MODE" "${CHROM_SIZE_INPUTS[@]}" <<'PY'
import sys
from pathlib import Path
import pysam
from nucleosuite.core.bam_headers import merge_bam_reference_headers_with_aliases
from nucleosuite.core.chrom_sizes import read_chrom_sizes_source, write_chrom_sizes_table
from nucleosuite.core.fragment_inputs import IntervalFragmentSource
from nucleosuite.core.regions import expand_contig_tokens

source_path, fasta_path, all_out, selected_out, tokens_text, input_mode, *inputs = sys.argv[1:]
if source_path:
    rows = list(read_chrom_sizes_source(source_path))
elif input_mode == "bam":
    handles = [pysam.AlignmentFile(path, "rb") for path in inputs]
    try:
        merged = merge_bam_reference_headers_with_aliases(handles)
        rows = list(zip(merged.references, merged.lengths))
    finally:
        for handle in handles:
            handle.close()
else:
    fasta = pysam.FastaFile(fasta_path)
    try:
        source = IntervalFragmentSource(inputs, fasta=fasta)
        try:
            rows = list(zip(source.references, source.lengths))
        finally:
            source.close()
    finally:
        fasta.close()

write_chrom_sizes_table(rows, all_out)
references = [name for name, _length in rows]
lengths = dict(rows)
selected = [
    value for value in expand_contig_tokens(tokens_text.split() or ["autosomes"], references)
    if ":" not in value
]
if not selected:
    raise SystemExit("No compatible contigs were selected for analysis")
write_chrom_sizes_table([(name, lengths[name]) for name in selected], selected_out)
PY
resolve_effective_blacklist() {
    if [[ "$NO_BLACKLIST" -eq 1 ]]; then
        BLACKLIST_BED=""
        echo "[INFO] Blacklist filtering disabled"
    elif [[ -n "$BLACKLIST_BED" ]]; then
        BLACKLIST_BED="$(realpath "$BLACKLIST_BED")"
        echo "[INFO] Using blacklist override: $BLACKLIST_BED"
    elif "$PYTHON_BIN" - "$CHROM_SIZES" <<'PY'
import sys
from nucleosuite.core.blacklist import is_hg19_reference
from nucleosuite.core.chrom_sizes import read_chrom_sizes_source
rows = list(read_chrom_sizes_source(sys.argv[1]))
raise SystemExit(0 if is_hg19_reference(dict(rows)) else 1)
PY
    then
        BLACKLIST_BED="$("$NUCLEOSUITE_BIN" resources path hg19-blacklist-v2)"
        BLACKLIST_BED="$(realpath "$BLACKLIST_BED")"
        echo "[INFO] Detected hg19/GRCh37 lengths; using bundled hg19 blacklist v2"
    else
        BLACKLIST_BED=""
        echo "[INFO] Reference is not confirmed hg19/GRCh37; no default blacklist applied"
    fi
    BLACKLIST_ARGS=()
    [[ -n "$BLACKLIST_BED" ]] && BLACKLIST_ARGS=(--blacklist-bed "$BLACKLIST_BED")
    return 0
}
resolve_effective_blacklist
{
    printf 'EFFECTIVE_BLACKLIST_BED\t%s\n' "$BLACKLIST_BED"
    printf 'BLACKLIST_MODE\t%s\n' "$([[ "$NO_BLACKLIST" -eq 1 ]] && echo disabled || { [[ -n "$BLACKLIST_BED" ]] && echo enabled || echo unavailable; })"
} >> "$PARAMETERS"
EFFECTIVE_IDENTITY_ARGS=(SELECTED_CHROM_SIZES "$CHROM_SIZES")
[[ -z "$BLACKLIST_BED" ]] || EFFECTIVE_IDENTITY_ARGS+=(EFFECTIVE_BLACKLIST "$BLACKLIST_BED")
append_input_identities "$PARAMETERS" "${EFFECTIVE_IDENTITY_ARGS[@]}"
PARAM_HASH="$("$PYTHON_BIN" - "$PARAMETERS" <<'PY'
import hashlib, sys
print(hashlib.sha256(open(sys.argv[1], 'rb').read()).hexdigest()[:16])
PY
)"

RANDOM_BED=""
ACTIVE_MAX_DUPLICATES="$MAX_DUPLICATES"
if [[ "$RUN_MODE" == "randomized" ]]; then
    if [[ "$RANDOMIZED_INPUT_READY" -eq 1 ]]; then
        [[ "$INPUT_MODE" == "fragments" && "${#FRAGMENTS[@]}" -eq 1 ]] || fatal "internal randomized-control input requires exactly one fragment BED"
        RANDOM_BED="${FRAGMENTS[0]}"
    else
        RANDOM_LOWER="$PNS_FRAG_LOWER"
        RANDOM_UPPER="$PNS_FRAG_UPPER"
        for length in "${EXACT_LENGTHS[@]}"; do
            (( length < RANDOM_LOWER )) && RANDOM_LOWER="$length"
            (( length > RANDOM_UPPER )) && RANDOM_UPPER="$length"
        done
        for spec in "${RANGE_SPECS[@]}"; do
            lo="${spec%%:*}"; hi="${spec##*:}"
            (( lo < RANDOM_LOWER )) && RANDOM_LOWER="$lo"
            (( hi > RANDOM_UPPER )) && RANDOM_UPPER="$hi"
        done
        (( FRAG_COUNT_MIN < RANDOM_LOWER )) && RANDOM_LOWER="$FRAG_COUNT_MIN"
        (( FRAG_COUNT_MAX > RANDOM_UPPER )) && RANDOM_UPPER="$FRAG_COUNT_MAX"
        (( HEATMAP_MIN_FRAG < RANDOM_LOWER )) && RANDOM_LOWER="$HEATMAP_MIN_FRAG"
        (( HEATMAP_MAX_FRAG > RANDOM_UPPER )) && RANDOM_UPPER="$HEATMAP_MAX_FRAG"
        RANDOM_BASE="$SETUP_DIR/$SAMPLE"
        RANDOM_BED="${RANDOM_BASE}.randomized.fragments.bed.gz"
        run_step "00_randomize_fragments" "$RANDOM_BED" "$NUCLEOSUITE_BIN" randomize-fragments \
            "${SOURCE_INPUT_ARGS[@]}" "${BLACKLIST_ARGS[@]}" --fasta "$FASTA" --chrom-sizes "$CHROM_SIZES" \
            -c "${CONTIGS[@]}" -o "$RANDOM_BASE" --frag-lower "$RANDOM_LOWER" --frag-upper "$RANDOM_UPPER" \
            --max-duplicates "$MAX_DUPLICATES" --dedup-scope "$DEDUP_SCOPE" --method dinucleotide \
            --seed "$RANDOMIZE_SEED" --search-window "$RANDOMIZE_SEARCH_WINDOW" --fallback "$RANDOMIZE_FALLBACK"
        RANDOM_MARKER="$DONE_DIR/${SUPPORT_PREFIX}${PARAM_HASH}_00_randomize_fragments.done"
        [[ -f "$RANDOM_MARKER" && -s "$RANDOM_BED" ]] || fatal "Randomized fragment generation failed; downstream analysis was not started"
    fi
    ANALYSIS_INPUT_ARGS=(--fragments "$RANDOM_BED")
    SOURCE_REFERENCE_ARGS=(--fasta "$FASTA")
    FRAGMENTS=("$RANDOM_BED")
    INPUT_MODE="fragments"
    # Source duplicate filtering was already applied before randomization.
    # Retain every materialized randomized row, including chance collisions.
    ACTIVE_MAX_DUPLICATES=0
fi

CTCF_FILTERED="$SETUP_DIR/${SUPPORT_PREFIX}ctcf_compatible.bed"
STATES_FILTERED="$SETUP_DIR/${SUPPORT_PREFIX}states_compatible.bed"
CTCF_EXPANDED="$SETUP_DIR/${SUPPORT_PREFIX}ctcf_compatible_flank${CTCF_FLANK}.bed"
GENE_SET_SUMMARY=""
GENE_STATE_INTERVAL=""
if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 0 ]]; then

# Reject cached CTCF files that retain motif strand in column 4 rather than
# the BED6 strand field in column 6.
if [[ -s "$CTCF_FILTERED" ]] && ! "$PYTHON_BIN" - "$CTCF_FILTERED" <<'PY'
import sys
valid = 0
with open(sys.argv[1], encoding="utf-8-sig") as handle:
    for raw in handle:
        text = raw.strip()
        if not text or text.startswith(("#", "track", "browser")):
            continue
        fields = text.split()
        if len(fields) < 6 or fields[5] not in {"+", "-"}:
            raise SystemExit(1)
        valid += 1
raise SystemExit(0 if valid else 1)
PY
then
  echo "[INFO] Rebuilding CTCF intervals as strand-aware BED6"
  rm -f "$CTCF_FILTERED"
  # Existing aggregate products were generated without valid motif orientation.
  # Remove them so resume mode recalculates every CTCF-centred profile.
  find "$AGG_DIR" -type f -delete
fi

run_step "00_prepare_regions" "$CTCF_FILTERED" "$PYTHON_BIN" - "$CHROM_SIZES" "$CTCF_BED" "$STATES_BED" "$CTCF_FILTERED" "$STATES_FILTERED" <<'PY'
import sys
from pathlib import Path
sizes_path,ctcf_path,states_path,ctcf_out,states_out=sys.argv[1:]
lengths={}
with open(sizes_path, encoding='utf-8-sig') as handle:
    for raw in handle:
        fields=raw.split()
        if len(fields) >= 2:
            lengths[fields[0]]=int(fields[1])
if not lengths:
    raise SystemExit(f'No chromosome sizes were found in {sizes_path}')
selected=list(lengths)
selected_set=set(selected)
def resolve(name):
    if name in selected_set: return name
    stripped=name[3:] if name.lower().startswith('chr') else name
    matches=[x for x in (stripped,'chr'+stripped) if x in selected_set]
    return matches[0] if len(matches)==1 else None
def compatible_rows(source):
    if not source:
        return
    with open(source,encoding='utf-8-sig') as inp:
        for raw in inp:
            text=raw.strip()
            if not text or text.startswith(('#','track','browser')): continue
            fields=text.split('\t') if '\t' in text else text.split()
            if len(fields)<3: continue
            chrom=resolve(fields[0])
            if chrom is None: continue
            try: start=max(0,int(fields[1])); end=min(lengths[chrom],int(fields[2]))
            except ValueError: continue
            if end<=start: continue
            yield chrom,start,end,fields

def filter_ctcf(source,dest):
    records=[]
    for chrom,start,end,fields in compatible_rows(source):
        strand=None
        for index in (5,3):
            if index < len(fields) and fields[index] in ('+','-'):
                strand=fields[index]
                break
        if strand is None:
            continue
        name = fields[3] if len(fields)>3 and fields[3] not in ('+','-','.') else f'CTCF_{len(records)+1}'
        score='0'
        for index in (4,8):
            if index < len(fields):
                try:
                    score=str(min(1000,max(0,int(round(float(fields[index]))))))
                    break
                except ValueError:
                    pass
        records.append((chrom,start,end,name,score,strand))
    if not records:
        raise SystemExit('No compatible stranded CTCF records remain')
    with open(dest,'w') as out:
        for record in records:
            out.write('\t'.join(map(str,record))+'\n')

def filter_bed(source,dest):
    if not source:
        Path(dest).write_text(''); return
    with open(dest,'w') as out:
        for chrom,start,end,fields in compatible_rows(source):
            fields[0]=chrom; fields[1]=str(start); fields[2]=str(end)
            out.write('\t'.join(fields)+'\n')

filter_ctcf(ctcf_path,ctcf_out)
filter_bed(states_path,states_out)
PY

run_step "00_expand_ctcf_regions" "$CTCF_EXPANDED" "$PYTHON_BIN" - \
  "$CTCF_FILTERED" "$CHROM_SIZES" "$CTCF_EXPANDED" "$CTCF_FLANK" <<'PY'
import sys
from pathlib import Path

source, sizes_path, output, flank_text = sys.argv[1:]
flank = int(flank_text)
sizes = {}
with open(sizes_path, encoding="utf-8") as handle:
    for raw in handle:
        chrom, length = raw.rstrip("\n").split("\t")[:2]
        sizes[chrom] = int(length)
rows = []
with open(source, encoding="utf-8") as handle:
    for raw in handle:
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 6 or fields[0] not in sizes:
            continue
        start, end = int(fields[1]), int(fields[2])
        center = (start + end) // 2
        left = max(0, center - flank)
        right = min(sizes[fields[0]], center + flank + 1)
        if right > left:
            rows.append([fields[0], str(left), str(right), *fields[3:6]])
Path(output).write_text("".join("\t".join(row) + "\n" for row in rows))
if not rows:
    raise SystemExit("No CTCF intervals could be expanded")
PY

if [[ -n "$GENES_BED" && -n "$STATES_BED" ]]; then
  GENE_SET_OUTPUT_PREFIX="gene_sets"
  GENE_SET_MEMBER_ARGS=()
  if [[ "$RUN_MODE" == "randomized" ]]; then
    GENE_SET_OUTPUT_PREFIX="${SAMPLE}_gene_sets"
    GENE_SET_MEMBER_ARGS=(--prefix-member-files)
  fi
  GENE_SET_SUMMARY="$GENE_SET_DIR/${GENE_SET_OUTPUT_PREFIX}_summary.tsv"
  # Rebuild any gene-set output whose summary does not reference the
  # strand-aware per-category TSS interval files required below.
  if [[ -s "$GENE_SET_SUMMARY" ]] && \
     ! head -n 1 "$GENE_SET_SUMMARY" | grep -q $'\tfinal_tss_interval'; then
    echo "[INFO] Rebuilding gene sets with strand-aware TSS intervals"
    rm -f "$GENE_SET_SUMMARY"
    rm -rf "$GENE_SET_DIR/final_tss" "$TSS_AGG_DIR"
    mkdir -p "$TSS_AGG_DIR"
  fi
  IFS=',' read -r -a VENN_SET_ARRAY <<< "$VENN_SETS"
  run_step "00_gene_sets" "$GENE_SET_SUMMARY" "$NUCLEOSUITE_BIN" gene-sets \
    "${BLACKLIST_ARGS[@]}" --genes-bed "$GENES_BED" --states-bed "$STATES_BED" --config "$GENE_SET_CONFIG" \
    --state-label-column "$STATES_LABEL_COLUMN" --chrom-sizes "$CHROM_SIZES" \
    --output-dir "$GENE_SET_DIR" --output-prefix "$GENE_SET_OUTPUT_PREFIX" --leftover-set-name leftover_genes \
    --venn-sets "${VENN_SET_ARRAY[@]}" "${GENE_SET_MEMBER_ARGS[@]}" --interval-format "$INTERVAL_FORMAT"
  GENE_STATE_INTERVAL="$GENE_SET_DIR/${GENE_SET_OUTPUT_PREFIX}_final_states.bed"
  [[ "$INTERVAL_FORMAT" == "bigbed" ]] && GENE_STATE_INTERVAL="$GENE_SET_DIR/${GENE_SET_OUTPUT_PREFIX}_final_states.bb"
fi
fi

PNS_BASE="$PNS_DIR/${SAMPLE}_PNS"
PNS_PREFIX="${PNS_BASE}_methodpns_mode${PNS_MODE_LENGTH}_lower${PNS_FRAG_LOWER}_upper${PNS_FRAG_UPPER}_smooth${PNS_SMOOTH_WINDOW}x${PNS_SMOOTH_ORDER}"
PNS_BW="${PNS_PREFIX}_pns.bw"
PNS_TRACK_LIST="pns,posPNS,coverage,dyad,fragment_ends,fragment_left_ends,fragment_right_ends,pns_peaks"
if [[ "$PNS_SMOOTH_WINDOW" -gt 0 ]]; then
  PNS_BW="${PNS_PREFIX}_pns_smoothed.bw"
  PNS_TRACK_LIST="pns_smoothed,${PNS_TRACK_LIST}"
fi
PNS_COVERAGE_BW="${PNS_PREFIX}_coverage.bw"
PNS_POS_BW="${PNS_PREFIX}_posPNS.bw"
PNS_COVERAGE_SCALED_BW="$SCALED_DIR/${SAMPLE}_coverage_mean_scaled.bw"
PNS_ANALYSIS_BW="$PNS_BW"


declare -A EXACT_DYADS EXACT_LEFT EXACT_RIGHT RANGE_DYADS RANGE_LEFT RANGE_RIGHT TYPE_TRACKS
WW_TYPE_LENGTH_TABLES=()
TRACK_SPEC_FILE="$COMBINED_TRACK_DIR/${SUPPORT_PREFIX}manifest.tsv"
printf 'fragment_range\toutput_prefix\ttracks\tbasic_scope\n' > "$TRACK_SPEC_FILE"
printf '%s\t%s\t%s\trange\n' "${PNS_FRAG_LOWER}-${PNS_FRAG_UPPER}" "$PNS_PREFIX" "$PNS_TRACK_LIST" >> "$TRACK_SPEC_FILE"

for length in "${EXACT_LENGTHS[@]}"; do
  folder="$EXACT_DIR/$length"; base="$folder/${SAMPLE}_${length}"
  prefix="${base}_dyads_lower${length}_upper${length}"
  EXACT_DYADS[$length]="${prefix}_dyad.bw"
  printf '%s\t%s\tdyad\trange\n' "$length" "$prefix" >> "$TRACK_SPEC_FILE"
  endfolder="$END_EXACT_DIR/$length"; endprefix="$endfolder/${SAMPLE}_${length}_fragment_ends_lower${length}_upper${length}"
  EXACT_LEFT[$length]="${endprefix}_fragment_left_ends.bw"; EXACT_RIGHT[$length]="${endprefix}_fragment_right_ends.bw"
  printf '%s\t%s\tfragment_left_ends,fragment_right_ends\trange\n' "$length" "$endprefix" >> "$TRACK_SPEC_FILE"
  dinuc_prefix="$SEQ_DIR/dinucleotide_profiles/exact/$length/${SAMPLE}_exact_${length}_dinuc_lower${length}_upper${length}"
  printf '%s\t%s\tdinuc_profile\trange\n' "$length" "$dinuc_prefix" >> "$TRACK_SPEC_FILE"
done

for spec in "${RANGE_SPECS[@]}"; do
  lo="${spec%%:*}"; hi="${spec##*:}"; label="${lo}_${hi}"; dir_label="${lo}-${hi}"
  prefix="$RANGE_DIR/$dir_label/${SAMPLE}_${label}_dyads_lower${lo}_upper${hi}"
  RANGE_DYADS[$label]="${prefix}_dyad.bw"
  printf '%s-%s\t%s\tdyad\trange\n' "$lo" "$hi" "$prefix" >> "$TRACK_SPEC_FILE"
  endprefix="$END_RANGE_DIR/$dir_label/${SAMPLE}_${label}_fragment_ends_lower${lo}_upper${hi}"
  RANGE_LEFT[$label]="${endprefix}_fragment_left_ends.bw"; RANGE_RIGHT[$label]="${endprefix}_fragment_right_ends.bw"
  printf '%s-%s\t%s\tfragment_left_ends,fragment_right_ends\trange\n' "$lo" "$hi" "$endprefix" >> "$TRACK_SPEC_FILE"
  dinuc_prefix="$SEQ_DIR/dinucleotide_profiles/ranges/$dir_label/${SAMPLE}_range_${label}_dinuc_lower${lo}_upper${hi}"
  printf '%s-%s\t%s\tdinuc_profile\trange\n' "$lo" "$hi" "$dinuc_prefix" >> "$TRACK_SPEC_FILE"
  ww_prefix="$SEQ_DIR/ww_types/ranges/$dir_label/${SAMPLE}_${label}_wwtypes_lower${lo}_upper${hi}"
  type_prefix="$SEQ_DIR/type_dyads/ranges/$dir_label/${SAMPLE}_${label}_wwtypes_lower${lo}_upper${hi}"
  printf '%s-%s\t%s\tww_types\trange\n' "$lo" "$hi" "$ww_prefix" >> "$TRACK_SPEC_FILE"
  printf '%s-%s\t%s\ttype_dyads\trange\n' "$lo" "$hi" "$type_prefix" >> "$TRACK_SPEC_FILE"
  WW_TYPE_LENGTH_TABLES+=("${ww_prefix}_ww_type_by_length.tsv")
  for type in 1 2 3 4; do TYPE_TRACKS["${label}_type${type}"]="${type_prefix}_type${type}_dyad.bw"; done
done

TRACK_REPORT="$COMBINED_TRACK_DIR/${SUPPORT_PREFIX}completion_report.tsv"
TRACK_ARGS=(
  "${ANALYSIS_INPUT_ARGS[@]}" "${BLACKLIST_ARGS[@]}" --fasta "$FASTA" --chrom-sizes "$CHROM_SIZES" -c "${CONTIGS[@]}"
  --output-dir "$COMBINED_TRACK_DIR" --spec-file "$TRACK_SPEC_FILE" --max-duplicates "$ACTIVE_MAX_DUPLICATES" --max-per-coordinate "$MAX_PER_COORDINATE"
  --dedup-scope "$DEDUP_SCOPE" --even-dyad "$EVEN_DYAD" --score-mode-length "$PNS_MODE_LENGTH" --bigbed-score-scale "$BIGBED_SCORE_SCALE"
  --score-smooth-window "$PNS_SMOOTH_WINDOW" --score-smooth-order "$PNS_SMOOTH_ORDER" --score-max-neg-run "$PNS_MAX_NEG_RUN"
  --interval-format "$INTERVAL_FORMAT" --output-format bigwig --report "$TRACK_REPORT"
)
run_step "01_combined_tracks" "$TRACK_REPORT" "$NUCLEOSUITE_BIN" tracks "${TRACK_ARGS[@]}"

PNS_NUC_RAW="${PNS_PREFIX}_nucleosome_regions.${INTERVAL_EXT}"
PNS_BRK_RAW="${PNS_PREFIX}_breakpoint_peaks.${INTERVAL_EXT}"
PNS_NUC="$PNS_NUC_RAW"
PNS_BRK="$PNS_BRK_RAW"

if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 0 ]]; then
  mkdir -p "$SCALED_DIR"
  run_step "01_scale_coverage" "$PNS_COVERAGE_SCALED_BW" "$NUCLEOSUITE_BIN" mean-scale \
    "$PNS_COVERAGE_BW" --scale 100 --output "$PNS_COVERAGE_SCALED_BW"
fi

RANGE_LABELS=()
for spec in "${RANGE_SPECS[@]}"; do RANGE_LABELS+=("${spec/:/_}"); done

if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 0 ]]; then

# Sequence profiles, WW/SS classifications and type-specific dyads were all
# generated by 01_combined_tracks. Summarise selected exact lengths from those tables.
EXACT_WW_TYPE_DIR="$SEQ_DIR/summaries"; mkdir -p "$EXACT_WW_TYPE_DIR"
EXACT_WW_TYPE_TSV="$EXACT_WW_TYPE_DIR/${SAMPLE}_exact_ww_type_by_length.tsv"
EXACT_WW_TYPE_PNG="$EXACT_WW_TYPE_DIR/${SAMPLE}_exact_ww_type_by_length_stacked.${PLOT_EXT}"
run_step "01_ww_type_exact_length_comparison" "$EXACT_WW_TYPE_PNG" "$PYTHON_BIN" - \
  "$EXACT_WW_TYPE_TSV" "$EXACT_WW_TYPE_PNG" "$(IFS=,; echo "${EXACT_LENGTHS[*]}")" \
  "${WW_TYPE_LENGTH_TABLES[@]}" <<'PYWW'
import sys
from nucleosuite.profile_plots import plot_ww_type_length_stacked
from nucleosuite.sequence.ww_types import write_selected_length_summary
output_tsv,output_png,selected_text,*input_tables=sys.argv[1:]
selected=[int(v) for v in selected_text.split(',') if v]
write_selected_length_summary(input_tables,selected,output_tsv)
plot_ww_type_length_stacked(output_tsv,output_png,title='WW/SS type frequencies for exact cfDNA fragment lengths')
PYWW

# 02_dac: DAC is calculated only from ranged dyad tracks.
run_dac_scopes() {
  local step_label="$1" track="$2" output_root="$3" output_label="$4"
  local genome_dir="$output_root/combined_chromosomes"
  mkdir -p "$genome_dir"
  queue_memory_step "02_dac_${step_label}_combined_chromosomes" "$genome_dir/${SAMPLE}_${output_label}_combined_chromosomes*_DAC_*.tsv" \
    "$NUCLEOSUITE_BIN" dac "${BLACKLIST_ARGS[@]}" --bigwig "$track" --chrom-sizes "$CHROM_SIZES" --scope combined_chromosomes \
    --window-size "$DAC_WINDOW_SIZE" --dmax "$DAC_DMAX" --algorithm "$DAC_ALGORITHM" \
    --out-prefix "${SAMPLE}_${output_label}_combined_chromosomes" --output-dir "$genome_dir" --progress-every 100
  if [[ -n "$GENE_STATE_INTERVAL" && -f "$GENE_STATE_INTERVAL" ]]; then
    local gene_dir="$output_root/gene_sets"
    mkdir -p "$gene_dir"
    queue_memory_step "02_dac_${step_label}_gene_sets" "$gene_dir/${SAMPLE}_${output_label}_gene_sets*_DAC_*.tsv" \
      "$NUCLEOSUITE_BIN" dac "${BLACKLIST_ARGS[@]}" --bigwig "$track" --regions-bed "$GENE_STATE_INTERVAL" \
      --state-column 4 --state-name gene_sets --dmax "$DAC_DMAX" --algorithm "$DAC_ALGORITHM" \
      --out-prefix "${SAMPLE}_${output_label}_gene_sets" --output-dir "$gene_dir" --progress-every 100
  fi
}
for spec in "${RANGE_SPECS[@]}"; do
  lo="${spec%%:*}"; hi="${spec##*:}"; label="${lo}_${hi}"; dir_label="${lo}-${hi}"
  run_dac_scopes "dyad_range_${label}" "${RANGE_DYADS[$label]}" \
    "$DAC_DIR/dyads/ranges/$dir_label" "dyad_${label}"
done
wait_queued_steps

DAC_VALIDATION="$DAC_DIR/${SUPPORT_PREFIX}DAC_COLUMN_VALIDATION.tsv"
run_step "02_verify_dac_columns" "$DAC_VALIDATION" "$PYTHON_BIN" - "$DAC_DIR" "$DAC_VALIDATION" <<'PYDAC'
import csv,sys
from pathlib import Path
root=Path(sys.argv[1]); output=Path(sys.argv[2]); rows=[]
required={'Distance','DAC Value','DAC Value Percent','Raw DAC Value','Opportunities'}
for path in sorted(root.rglob('*DAC*.tsv')):
    if path.name.endswith('_summary.tsv') or path == output:
        continue
    with path.open(encoding='utf-8-sig') as handle:
        header=set(next(csv.reader(handle,delimiter='\t'),[]))
    missing=sorted(required-header)
    rows.append((str(path), 'PASS' if not missing else 'FAIL', ','.join(missing)))
if not rows:
    raise SystemExit('No DAC result TSVs found')
output.write_text('path\tstatus\tmissing_columns\n'+''.join(f'{p}\t{s}\t{m}\n' for p,s,m in rows))
if any(s=='FAIL' for _,s,_ in rows):
    raise SystemExit('One or more DAC files lack required columns')
PYDAC

# 04_nrl: long-, short-, and nucleosome-scale periodicity from ranged-dyad DAC curves.
run_nrl_analysis() {
  local input_tsv="$1" output_dir="$2" analysis_label="$3" min_distance="$4" max_distance="$5" peak_resolution="$6" skip_first="${7:-0}"
  local input_stem safe_stem output_prefix parameter_prefix step_name
  mkdir -p "$output_dir"
  input_stem="$(basename "$input_tsv" .tsv)"; safe_stem="${input_stem//[^A-Za-z0-9._-]/_}"
  output_prefix="$output_dir/${input_stem}_${analysis_label}"
  parameter_prefix="${output_prefix}_peakres${peak_resolution}_min${min_distance}_max${max_distance}_skipfirst${skip_first}"
  step_name="04_nrl_${safe_stem}_${analysis_label}"
  queue_step "$step_name" "${parameter_prefix}_regression.tsv" "$NUCLEOSUITE_BIN" nrl "$input_tsv" \
    --min-distance "$min_distance" --max-distance "$max_distance" --peak-resolution "$peak_resolution" \
    --skip-first-peaks "$skip_first" --output-prefix "$output_prefix" --title "$input_stem"
}
run_nrl_tree() {
  local source_root="$1" target_root="$2" range_hi="$3"
  local input_tsv relative parent outdir intermediate_min
  intermediate_min=$((range_hi + 1))
  while IFS= read -r -d '' input_tsv; do
    [[ "$input_tsv" == *_summary.tsv ]] && continue
    [[ "$(basename "$input_tsv")" == *DAC_COLUMN_VALIDATION.tsv ]] && continue
    relative="${input_tsv#"$source_root"/}"; parent="$(dirname "$relative")"; outdir="$target_root/$parent"
    run_nrl_analysis "$input_tsv" "$outdir" "nrl_${NRL_MIN_DISTANCE}_${NRL_MAX_DISTANCE}" \
      "$NRL_MIN_DISTANCE" "$NRL_MAX_DISTANCE" "$NRL_PEAK_RESOLUTION" 1
    run_nrl_analysis "$input_tsv" "$outdir" "periodicity_${SHORT_PERIODICITY_MIN}_${SHORT_PERIODICITY_MAX}" \
      "$SHORT_PERIODICITY_MIN" "$SHORT_PERIODICITY_MAX" 1 0
    run_nrl_analysis "$input_tsv" "$outdir" "periodicity_${intermediate_min}_${INTERMEDIATE_PERIODICITY_MAX}" \
      "$intermediate_min" "$INTERMEDIATE_PERIODICITY_MAX" "$INTERMEDIATE_PERIODICITY_RESOLUTION" 0
  done < <(find "$source_root" -type f -name "*DAC*.tsv" -print0 | sort -z)
}
if [[ "$SKIP_NRL" -eq 0 ]]; then
  for spec in "${RANGE_SPECS[@]}"; do
    lo="${spec%%:*}"; hi="${spec##*:}"; dir_label="${lo}-${hi}"
    run_nrl_tree "$DAC_DIR/dyads/ranges/$dir_label" "$NRL_DIR/from_dac/dyads/ranges/$dir_label" "$hi"
  done
fi
wait_queued_steps

# 05_ctcf_aggregation: every signal is routed to its own input-track directory.
aggregate_track() {
  local step="$1" label="$2" track="$3" output_dir="$4" ylabel="$5"
  local -a sparse_args=(--zero-thresh 0 --nan-to-zero)
  [[ "$label" == dyad_* || "$label" == type* ]] && sparse_args+=(--max-score 1000000000000)
  mkdir -p "$output_dir"
  local prefix="${SAMPLE}_CTCF_${label}"
  queue_step "$step" "$output_dir/${prefix}_win*_heatmap.${PLOT_EXT}" "$NUCLEOSUITE_BIN" aggregate \
    "${BLACKLIST_ARGS[@]}" --bigwig "$track" --region-bed "$CTCF_FILTERED" --output-dir "$output_dir" --output-prefix "$prefix" \
    --window-half "$AGGREGATE_WINDOW_HALF" --strand-col 6 --missing-strand error "${sparse_args[@]}" \
    --mean-ylabel "$ylabel" --colorbar-label "$ylabel"
}
aggregate_track "05_ctcf_pns" pns "$PNS_ANALYSIS_BW" "$AGG_DIR/pns" "Mean PNS"
for length in "${EXACT_LENGTHS[@]}"; do
  aggregate_track "05_ctcf_dyad_exact_${length}" "dyad_exact_${length}" "${EXACT_DYADS[$length]}" \
    "$AGG_DIR/dyads/exact/$length" "Mean ${length} bp dyad signal"
done
for spec in "${RANGE_SPECS[@]}"; do
  lo="${spec%%:*}"; hi="${spec##*:}"; label="${lo}_${hi}"; dir_label="${lo}-${hi}"
  aggregate_track "05_ctcf_dyad_range_${label}" "dyad_range_${label}" "${RANGE_DYADS[$label]}" \
    "$AGG_DIR/dyads/ranges/$dir_label" "Mean ${dir_label} bp dyad signal"
  for type in 1 2 3 4; do
    aggregate_track "05_ctcf_type_dyad_${label}_type${type}" "type${type}_${label}" \
      "${TYPE_TRACKS[${label}_type${type}]}" "$AGG_DIR/type_dyads/ranges/$dir_label/type${type}" \
      "Mean type ${type} dyad signal"
  done
done
wait_queued_steps

# 06_tss_aggregation: organised first by input track, then by gene set.
tss_aggregate_track() {
  local label="$1" track="$2" track_dir="$3" ylabel="$4"
  local summary="$GENE_SET_SUMMARY"
  [[ -s "$summary" ]] || return 0
  mkdir -p "$track_dir"
  local set_name tss_interval prefix profile mean_png set_dir
  local -a profile_specs=() sparse_args=(--zero-thresh 0 --nan-to-zero)
  [[ "$label" == dyad_* || "$label" == type* ]] && sparse_args+=(--max-score 1000000000000)
  while IFS=$'\t' read -r set_name final_count tss_interval; do
    [[ -n "$set_name" ]] || continue
    [[ "$final_count" =~ ^[0-9]+$ && "$final_count" -gt 0 ]] || continue
    [[ -n "$tss_interval" && -e "$tss_interval" ]] || continue
    set_dir="$track_dir/$set_name"; mkdir -p "$set_dir"
    prefix="${SAMPLE}_${label}_TSS_${set_name}"
    profile="$set_dir/${prefix}_aggregate_all.tsv"; mean_png="$set_dir/${prefix}_mean.${PLOT_EXT}"
    queue_step "06_tss_${label}_${set_name}" "$profile" "$NUCLEOSUITE_BIN" aggregate \
      "${BLACKLIST_ARGS[@]}" --bigwig "$track" --region-bed "$tss_interval" --output-dir "$set_dir" --output-prefix "$prefix" \
      --aggregate-output "$profile" --mean-plot-output "$mean_png" --window-half "$AGGREGATE_WINDOW_HALF" --strand-col 6 \
      --missing-strand error "${sparse_args[@]}" --mean-ylabel "$ylabel" --colorbar-label "$ylabel"
    profile_specs+=("${set_name}=${profile}")
  done < <("$PYTHON_BIN" - "$summary" <<'PYTSS'
import csv,sys
with open(sys.argv[1],encoding='utf-8') as handle:
    for row in csv.DictReader(handle,delimiter='\t'):
        print(row.get('set_name',''),row.get('final_gene_count','0'),row.get('final_tss_interval',''),sep='\t')
PYTSS
  )
  wait_queued_steps
  local -a available_profile_specs=()
  local profile_spec profile_path
  for profile_spec in "${profile_specs[@]}"; do
    profile_path="${profile_spec#*=}"
    [[ -s "$profile_path" ]] && available_profile_specs+=("$profile_spec")
  done
  profile_specs=("${available_profile_specs[@]}")
  if [[ ${#profile_specs[@]} -gt 0 ]]; then
    local combined_dir="$track_dir/combined_plots"; mkdir -p "$combined_dir"
    local combined_prefix="${SAMPLE}_${label}_TSS_gene_sets_combined"
    run_step "06_tss_overlay_${label}" "$combined_dir/${combined_prefix}.${PLOT_EXT}" "$PYTHON_BIN" - \
      "$combined_dir/${combined_prefix}.tsv" "$combined_dir/${combined_prefix}.${PLOT_EXT}" \
      "${SAMPLE}: ${label} at gene-set TSS" "$ylabel" "${profile_specs[@]}" <<'PYTSSPLOT'
import sys
from nucleosuite.profile_plots import plot_profile_overlay
output_tsv,output_png,title,ylabel,*specs=sys.argv[1:]
plot_profile_overlay([(s.split('=',1)[0],s.split('=',1)[1]) for s in specs],output_tsv,output_png,
                     xlabel='Position relative to TSS (bp)',ylabel=ylabel,title=title)
PYTSSPLOT
  fi
}
tss_aggregate_track pns "$PNS_ANALYSIS_BW" "$TSS_AGG_DIR/pns" "Mean PNS"
for length in "${EXACT_LENGTHS[@]}"; do
  tss_aggregate_track "dyad_exact_${length}" "${EXACT_DYADS[$length]}" "$TSS_AGG_DIR/dyads/exact/$length" "Mean ${length} bp dyad signal"
done
for spec in "${RANGE_SPECS[@]}"; do
  lo="${spec%%:*}"; hi="${spec##*:}"; label="${lo}_${hi}"; dir_label="${lo}-${hi}"
  tss_aggregate_track "dyad_range_${label}" "${RANGE_DYADS[$label]}" "$TSS_AGG_DIR/dyads/ranges/$dir_label" "Mean ${dir_label} bp dyad signal"
  for type in 1 2 3 4; do
    tss_aggregate_track "type${type}_${label}" "${TYPE_TRACKS[${label}_type${type}]}" \
      "$TSS_AGG_DIR/type_dyads/ranges/$dir_label/type${type}" "Mean type ${type} dyad signal"
  done
done


# 06_tss_expression_quintiles: TSS profiles stratified by tissue nTPM.
if [[ "$SKIP_TSS_EXPRESSION_QUINTILES" -eq 0 ]]; then
  [[ -n "$GENES_BED" ]] || fatal "TSS expression-quintile analysis requires --genes-bed or a resource set that supplies genes"
  TSS_TISSUE_KEY="${TSS_EXPRESSION_TISSUE// /_}"
  TSS_TISSUE_KEY="${TSS_TISSUE_KEY//[^A-Za-z0-9._-]/_}"
  run_tss_expression_quintiles() {
    local label="$1" signal="$2" output_dir="$3"
    mkdir -p "$output_dir"
    local prefix="$output_dir/${SAMPLE}_${label}_TSS_${TSS_TISSUE_KEY}_expression_quintiles"
    queue_step "06_tss_expression_quintiles_${label,,}" "${prefix}_window*_tss_expression_quintiles_metadata.tsv" \
      "$NUCLEOSUITE_BIN" tss-expression-quintiles --signal "$signal" --sample "$SAMPLE" \
      "${BLACKLIST_ARGS[@]}" --signal-label "$label" --expression "$TSS_EXPRESSION_RESOURCE" --tissue "$TSS_EXPRESSION_TISSUE" \
      --genes-bed "$GENES_BED" --window "$TSS_EXPRESSION_WINDOW" --output-prefix "$prefix"
  }
  run_tss_expression_quintiles PNS "$PNS_ANALYSIS_BW" "$TSS_EXPRESSION_DIR/$TSS_TISSUE_KEY/pns"
fi
wait_queued_steps

# 07_distances: adjacent spacing plus 1-7 order NRL regression from PNS nucleosome calls.
whole_dir="$DIST_DIR/pns_peaks/combined_chromosomes"; mkdir -p "$whole_dir"
adjacent_prefix="$whole_dir/${SAMPLE}_PNS_peak_distances_adjacent"
queue_memory_step "07_distances_pns_adjacent" "${adjacent_prefix}*.tsv" "$NUCLEOSUITE_BIN" distances "$PNS_NUC" \
  "${BLACKLIST_ARGS[@]}" --position-column 7 --score-column 5 --score-percentile 0 --min-distance 1 --max-distance "$DISTANCE_ADJACENT_MAX" \
  --max-order 1 --scope combined_chromosomes --write-filtered-bed --interval-format "$INTERVAL_FORMAT" \
  --interval-chrom-sizes "$CHROM_SIZES" --output-prefix "$adjacent_prefix"
long_prefix="$whole_dir/${SAMPLE}_PNS_peak_distances_orders1-${DISTANCE_LONG_MAX_ORDER}"
queue_memory_step "07_distances_pns_nrl" "${long_prefix}*.tsv" "$NUCLEOSUITE_BIN" distances "$PNS_NUC" \
  "${BLACKLIST_ARGS[@]}" --position-column 7 --score-column 5 --score-percentile 0 --min-distance 1 --max-distance "$DISTANCE_LONG_MAX" \
  --max-order "$DISTANCE_LONG_MAX_ORDER" --scope combined_chromosomes --regression-scope combined \
  --write-filtered-bed --interval-format "$INTERVAL_FORMAT" --interval-chrom-sizes "$CHROM_SIZES" --output-prefix "$long_prefix"
if [[ -n "$STATES_BED" ]]; then
  state_dir="$DIST_DIR/pns_peaks/chromhmm_states"; mkdir -p "$state_dir"
  state_prefix="$state_dir/${SAMPLE}_PNS_ChromHMM_peak_distances"
  queue_memory_step "07_state_distances_pns" "${state_prefix}_scorepct0_state_relative_percent.${PLOT_EXT}" \
    "$NUCLEOSUITE_BIN" distances "$PNS_NUC" --position-column 7 --score-column 5 --score-percentile 0 \
    "${BLACKLIST_ARGS[@]}" --min-distance 1 --max-distance "$STATE_DISTANCE_MAX" --max-order 1 --scope combined_chromosomes \
    --state-bed "$STATES_FILTERED" --state-label-column "$STATES_LABEL_COLUMN" --state-color-column 9 \
    --state-overlay-plot --state-overlay-smooth-window "$STATE_DISTANCE_SMOOTH_WINDOW" \
    --state-overlay-smooth-polyorder "$STATE_DISTANCE_SMOOTH_ORDER" \
    --state-overlay-title "${SAMPLE} PNS: adjacent peak distances by ChromHMM state" --output-prefix "$state_prefix"
fi
wait_queued_steps

# 08_region_extract: CTCF-centred regional tables organised by signal.
if [[ "$SKIP_REGION_EXTRACT" -eq 0 ]]; then
  pns_dir="$REGION_DIR/ctcf/pns"; mkdir -p "$pns_dir"
  PNS_REGION_PREFIX="$pns_dir/${SAMPLE}_CTCF_PNS"
  queue_step "08_region_extract_pns" "${PNS_REGION_PREFIX}_pns_signal.tsv" "$NUCLEOSUITE_BIN" region-extract \
    "${BLACKLIST_ARGS[@]}" --bed "$CTCF_EXPANDED" --coverage-bw "$PNS_COVERAGE_SCALED_BW" --score-bw "$PNS_ANALYSIS_BW" \
    --nucleosome-peaks "$PNS_NUC" --breakpoint-peaks "$PNS_BRK" --peak-flank-bp "$REGION_PEAK_FLANK" \
    --peak-center-column 7 --peak-score-column 5 --out-prefix "$PNS_REGION_PREFIX" \
    --chrom-mode auto --missing-chrom error --progress-every 100 --overwrite
fi
wait_queued_steps

# 09_fragment_lengths: organised by the region set counted.
WHOLE_FRAG_DIR="$FRAG_DIR/combined_chromosomes"; mkdir -p "$WHOLE_FRAG_DIR"
WHOLE_FRAG_TSV="$WHOLE_FRAG_DIR/${SAMPLE}_combined_chromosomes_fragment_lengths.tsv"
queue_step "09_fragment_lengths_combined_chromosomes" "$WHOLE_FRAG_TSV" "$NUCLEOSUITE_BIN" fragment-lengths "${ANALYSIS_INPUT_ARGS[@]}" "${BLACKLIST_ARGS[@]}" \
  --contigs "${CONTIGS[@]}" --min-length "$FRAG_COUNT_MIN" --max-length "$FRAG_COUNT_MAX" --output "$WHOLE_FRAG_TSV" \
  --plot "$WHOLE_FRAG_DIR/${SAMPLE}_combined_chromosomes_fragment_lengths.${PLOT_EXT}" --plot-min "$FRAG_PLOT_MIN" --plot-max "$FRAG_PLOT_MAX"
STATE_FRAG_TSV=""
if [[ -n "$STATES_BED" ]]; then
  state_frag_dir="$FRAG_DIR/chromhmm_states"; mkdir -p "$state_frag_dir"
  STATE_FRAG_TSV="$state_frag_dir/${SAMPLE}_states_fragment_lengths.tsv"
  queue_step "09_fragment_lengths_states" "$STATE_FRAG_TSV" "$NUCLEOSUITE_BIN" fragment-lengths "${ANALYSIS_INPUT_ARGS[@]}" "${BLACKLIST_ARGS[@]}" \
    --bed "$STATES_FILTERED" --bed-label-column "$STATES_LABEL_COLUMN" --overlap-policy all --contigs "${CONTIGS[@]}" \
    --min-length "$FRAG_COUNT_MIN" --max-length "$FRAG_COUNT_MAX" --output "$STATE_FRAG_TSV" --separate-files \
    --output-dir "$state_frag_dir/state_profiles" --plot "$state_frag_dir/${SAMPLE}_states_fragment_lengths.${PLOT_EXT}" \
    --plot-min "$FRAG_PLOT_MIN" --plot-max "$FRAG_PLOT_MAX"
fi
wait_queued_steps

# 10_fragment_heatmaps: combined visualisation of the available length profiles.
if [[ "$SKIP_FRAGMENT_HEATMAP" -eq 0 ]]; then
  missing_heatmap_dependencies=()
  [[ -s "$WHOLE_FRAG_TSV" ]] || missing_heatmap_dependencies+=("$WHOLE_FRAG_TSV")
  if [[ -n "$STATE_FRAG_TSV" && ! -s "$STATE_FRAG_TSV" ]]; then
    missing_heatmap_dependencies+=("$STATE_FRAG_TSV")
  fi
  if [[ "${#missing_heatmap_dependencies[@]}" -gt 0 ]]; then
    echo "[SKIP] 10_fragment_heatmap: missing fragment-length dependency: ${missing_heatmap_dependencies[*]}"
    SKIP_COUNT=$((SKIP_COUNT + 1))
  else
    heatmap_dir="$HEATMAP_DIR/combined"; mkdir -p "$heatmap_dir"
    args=(--input "Combined_chromosomes=$WHOLE_FRAG_TSV")
    [[ -z "$STATE_FRAG_TSV" ]] || args+=(--input "$STATE_FRAG_TSV")
    run_step "10_fragment_heatmap" "$heatmap_dir/${SAMPLE}_fragment_lengths_fragmin*_heatmap.${PLOT_EXT}" \
      "$NUCLEOSUITE_BIN" fragment-heatmap "${args[@]}" --out-prefix "$heatmap_dir/${SAMPLE}_fragment_lengths" \
      --min-frag "$HEATMAP_MIN_FRAG" --max-frag "$HEATMAP_MAX_FRAG" --normalization "$HEATMAP_NORMALIZATION" \
      --title "${SAMPLE}: cfDNA fragment-length profiles"
  fi
fi

# 11_gene_expression: grouped by the signal used for the analysis.
if [[ -n "$EXPRESSION" && "$SKIP_GENE_EXPRESSION" -eq 0 ]]; then
  FOCUS_PROFILE_ARGS=(); for profile in "${EXPRESSION_FOCUS_PROFILES[@]}"; do FOCUS_PROFILE_ARGS+=(--focus-profile "$profile"); done
  run_gene_expression() {
    local label="$1" signal_type="$2" peaks="$3" signal="$4"; local output_dir="$GENE_EXPRESSION_DIR/${label,,}"
    mkdir -p "$output_dir"; local prefix="$output_dir/${SAMPLE}_${label}_gene_expression"
    queue_memory_step "11_gene_expression_${label,,}" "${prefix}_analysis*_metadata.tsv" "$NUCLEOSUITE_BIN" gene-expression \
      "${BLACKLIST_ARGS[@]}" --expression "$EXPRESSION" --genes-bed "$GENES_BED" --peaks "${SAMPLE}=${peaks}" --signal "${SAMPLE}=${signal}" \
      --signal-type "$signal_type" --analysis all --output-prefix "$prefix" --expression-gene-column "$EXPRESSION_GENE_COLUMN" \
      --expression-name-column "$EXPRESSION_NAME_COLUMN" --expression-profile-column "$EXPRESSION_PROFILE_COLUMN" \
      --expression-value-column "$EXPRESSION_VALUE_COLUMN" --fft-window "$GENE_FFT_WINDOW" \
      --fft-period-min "$GENE_FFT_PERIOD_MIN" --fft-period-max "$GENE_FFT_PERIOD_MAX" \
      --fft-ranking-periods "$GENE_FFT_RANKING_PERIODS" "${FOCUS_PROFILE_ARGS[@]}"
  }
  run_gene_expression PNS pns "$PNS_NUC" "$PNS_ANALYSIS_BW"
fi
wait_queued_steps

# 12_positive_runs: the active suite input (observed or randomized-only).
if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 0 && "$SKIP_POSITIVE_RUNS" -eq 0 ]]; then
  pns_positive_dir="$POSITIVE_RUNS_DIR/pns"; mkdir -p "$pns_positive_dir"
  PNS_POSITIVE_PREFIX="$pns_positive_dir/${SAMPLE}_PNS_positive_runs"
  queue_step "12_positive_runs_pns" "${PNS_POSITIVE_PREFIX}_threshold*_summary.tsv" "$NUCLEOSUITE_BIN" positive-runs \
    "${BLACKLIST_ARGS[@]}" --bigwig "$PNS_BW" --output-prefix "$PNS_POSITIVE_PREFIX" --contigs "${CONTIGS[@]}" \
    --threshold "$POSITIVE_RUNS_THRESHOLD" --chunk-size "$POSITIVE_RUNS_CHUNK_SIZE" \
    --min-run-length "$POSITIVE_RUNS_MIN_LENGTH" --max-run-length "$POSITIVE_RUNS_MAX_LENGTH" \
    --plot-x-max "$POSITIVE_RUNS_PLOT_X_MAX" --normalization "$POSITIVE_RUNS_NORMALIZATION" \
    --title "${SAMPLE}: PNS positive run lengths"
fi
wait_queued_steps

# 13_peak_analysis: active-mode PNS peak-score distributions.
run_peak_score_frequency() {
  local step="$1" label="$2" peaks="$3" title="$4" output_root="$5"
  local output_dir="$output_root/$label"; mkdir -p "$output_dir"
  local prefix output
    prefix="$output_dir/${SAMPLE}_${label}"
    output="${prefix}_bins*_score_frequency.tsv"
  local -a peak_args=(--peaks "${RUN_MODE}=$peaks")
  queue_memory_step "$step" "$output" "$NUCLEOSUITE_BIN" peak-score-frequency "${peak_args[@]}" \
    "${BLACKLIST_ARGS[@]}" --output-prefix "$prefix" --score-column 5 --score-scale 1 --integer-bins \
    --normalization "$PEAK_SCORE_NORMALIZATION" --title "$title"
}
if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 0 && "$SKIP_PEAK_SCORE_FREQUENCY" -eq 0 ]]; then
  run_peak_score_frequency "13_peak_scores_pns_nucleosome" pns_nucleosome "$PNS_NUC" \
    "${SAMPLE}: PNS nucleosome-region scores" "$PEAK_ANALYSIS_DIR/pns/score_frequencies"
  run_peak_score_frequency "13_peak_scores_pns_breakpoint" pns_breakpoint "$PNS_BRK" \
    "${SAMPLE}: PNS breakpoint-peak scores" "$PEAK_ANALYSIS_DIR/pns/score_frequencies"
fi
wait_queued_steps

fi  # end analytical stages skipped by combine-prerequisites-only

# Guarantee that no asynchronous analysis remains before writing the final report.
wait_queued_steps

FAILED_STEPS_TSV="$OUTDIR/${SUPPORT_PREFIX}failed_steps.tsv"
if [[ "$FAIL_COUNT" -gt 0 ]]; then
  {
    echo -e "step\tstatus\treason\tlog"
    for index in "${!FAILED_STEPS[@]}"; do
      printf '%s\t%s\t%s\t%s\n' \
        "${FAILED_STEPS[$index]}" "${FAILED_STEP_STATUS[$index]}" \
        "${FAILED_STEP_REASON[$index]}" "${FAILED_STEP_LOG[$index]}"
    done
  } > "$FAILED_STEPS_TSV"
else
  rm -f "$FAILED_STEPS_TSV"
fi

REPORT="$OUTDIR/${SUPPORT_PREFIX}NUCLEOSUITE_CFDNA_SUITE_REPORT.tsv"
{
  echo -e "metric\tvalue"
  echo -e "sample\t$SAMPLE"
  echo -e "input_mode\t$INPUT_MODE"
  echo -e "pns_fragment_range\t${PNS_FRAG_LOWER}-${PNS_FRAG_UPPER}"
  echo -e "pns_mode_length\t$PNS_MODE_LENGTH"
  echo -e "bigbed_score_scale\t$BIGBED_SCORE_SCALE"
  echo -e "exact_lengths\t$(IFS=,; echo "${EXACT_LENGTHS[*]}")"
  echo -e "range_lengths\t$(IFS=,; echo "${RANGE_LABELS[*]}")"
  echo -e "gene_sets\tactive_genes,weak_genes,repressed_genes,leftover_genes"
  echo -e "expression_table\t${EXPRESSION:-}"
  echo -e "gene_expression_signals\tpns"
  echo -e "run_mode\t$RUN_MODE"
  echo -e "parameter_hash\t$PARAM_HASH"
  echo -e "blacklist_bed\t$BLACKLIST_BED"
  if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 1 ]]; then echo -e "execution_scope\tcombine_prerequisites_only"; else echo -e "execution_scope\tcombined_chromosomes_analysis"; fi
  echo -e "passed_steps\t$PASS_COUNT"
  echo -e "failed_steps\t$FAIL_COUNT"
  echo -e "skipped_completed_steps\t$SKIP_COUNT"
  if [[ "$FAIL_COUNT" -gt 0 ]]; then
    printf 'failed_step_names\t%s\n' "$(IFS=,; echo "${FAILED_STEPS[*]}")"
    echo -e "failed_steps_tsv\t$FAILED_STEPS_TSV"
  else
    echo -e "failed_step_names\t"
    echo -e "failed_steps_tsv\t"
  fi
} > "$REPORT"
MANIFEST="$OUTDIR/${SUPPORT_PREFIX}run_manifest.json"
"$PYTHON_BIN" - "$OUTDIR" "$PARAMETERS" "$REPORT" "$MANIFEST" "$SAMPLE" "$RUN_MODE" "$PARAM_HASH" "$BLACKLIST_BED" "$FAIL_COUNT" "$RANDOM_BED" <<'PYMANIFEST'
import csv, json, os, sys, tempfile
from pathlib import Path
from nucleosuite import __version__
root, parameters_path, report_path, output_path = map(Path, sys.argv[1:5])
sample, run_mode, parameter_hash, blacklist, failure_count, randomized_bed = sys.argv[5:]
parameters = {}
with parameters_path.open(encoding="utf-8") as handle:
    for row in csv.reader(handle, delimiter="\t"):
        if len(row) >= 2 and row[0] != "parameter":
            parameters[row[0]] = row[1]

def indexed_values(prefix):
    rows = []
    for key, value in parameters.items():
        if key.startswith(prefix + "_") and key[len(prefix) + 1:].isdigit():
            rows.append((int(key[len(prefix) + 1:]), value))
    return [value for _index, value in sorted(rows)]

input_identities = {}
for key, value in parameters.items():
    if key.startswith("INPUT_IDENTITY_"):
        try:
            input_identities[key.removeprefix("INPUT_IDENTITY_")] = json.loads(value)
        except json.JSONDecodeError:
            input_identities[key.removeprefix("INPUT_IDENTITY_")] = value
parameter_bams = indexed_values("BAM")
parameter_fragments = indexed_values("FRAGMENT")
active_bams = [] if randomized_bed else parameter_bams
active_fragments = (
    [str(Path(randomized_bed).resolve())]
    if randomized_bed
    else parameter_fragments
)
outputs = []
for path in sorted(root.rglob("*")):
    if not path.is_file() or ".done" in path.parts or path == output_path:
        continue
    outputs.append({"path": str(path.relative_to(root)), "size_bytes": path.stat().st_size})
payload = {
    "schema_version": 1, "nucleosuite_version": __version__, "suite": "cfdna",
    "sample": sample, "run_mode": run_mode, "parameter_hash": parameter_hash,
    "blacklist_bed": blacklist or None, "success": int(failure_count) == 0,
    "parameters_file": str(parameters_path.resolve()), "report_file": str(report_path.resolve()),
    "source_inputs": {
        "bam": indexed_values("SOURCE_BAM"),
        "fragments": indexed_values("SOURCE_FRAGMENT"),
    },
    "active_inputs": {
        "bam": active_bams,
        "fragments": active_fragments,
    },
    "input_identities": input_identities,
    "parameters": parameters, "outputs": outputs,
}
output_path.parent.mkdir(parents=True, exist_ok=True)
fd, temporary_name = tempfile.mkstemp(prefix=f".{output_path.name}.", suffix=".partial", dir=output_path.parent)
os.close(fd)
temporary = Path(temporary_name)
try:
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, output_path)
finally:
    temporary.unlink(missing_ok=True)
PYMANIFEST
echo "NucleoSuite cfDNA suite complete"
echo "Passed: $PASS_COUNT  Failed: $FAIL_COUNT  Skipped: $SKIP_COUNT"
echo "Report: $REPORT"
echo "Logs: $LOG_DIR"
if [[ "$FAIL_COUNT" -gt 0 ]]; then
  echo "Failed steps: $FAILED_STEPS_TSV" >&2
  printf '  - %s\n' "${FAILED_STEPS[@]}" >&2
  exit 1
fi
exit 0
