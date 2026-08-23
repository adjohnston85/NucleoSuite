#!/usr/bin/env bash
# Run a configurable, comprehensive NucleoSuite workflow for paired-end MNase BAM data.

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

PNS_FRAG_LOWER=120
PNS_FRAG_UPPER=180
PNS_MODE_LENGTH=147
BIGBED_SCORE_SCALE=1000
FINE_FRAG_LOWER=146
FINE_FRAG_UPPER=148
EXACT_SIZE=147
DINUC_EXACT_A=145
DINUC_EXACT_B=147
MAX_DUPLICATES=1
MAX_PER_COORDINATE=0
DEDUP_SCOPE="all_bams"
EVEN_DYAD="split"

PNS_SMOOTH_WINDOW=0
PNS_SMOOTH_ORDER=2
PNS_MAX_NEG_RUN=0
PEAK_SMOOTH_WINDOW=0
PEAK_SMOOTH_ORDER=2

CTCF_FLANK=2000
AGGREGATE_WINDOW_HALF=2500
REGION_PEAK_FLANK=2000
STATES_LABEL_COLUMN=4

DAC_DMAX=2000
DAC_WINDOW_SIZE=100000
DAC_ALGORITHM="auto"

NRL_MIN_DISTANCE=1
NRL_MAX_DISTANCE=1500
NRL_PEAK_RESOLUTION=160
DISTANCE_X_MAJOR_TICK=""
DISTANCE_X_MINOR_TICK=""
SHORT_PERIODICITY_MIN=1
SHORT_PERIODICITY_MAX=144
INTERMEDIATE_PERIODICITY_MIN=150
INTERMEDIATE_PERIODICITY_MAX=220
INTERMEDIATE_PERIODICITY_RESOLUTION=8

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

POSITIVE_RUNS_THRESHOLD=0
POSITIVE_RUNS_CHUNK_SIZE=1000000
POSITIVE_RUNS_MIN_LENGTH=1
POSITIVE_RUNS_MAX_LENGTH=0
POSITIVE_RUNS_PLOT_X_MAX=550
POSITIVE_RUNS_NORMALIZATION="count"

PEAK_SCORE_NORMALIZATION="count"

GENE_FFT_WINDOW=10000
GENE_FFT_PERIOD_MIN=120
GENE_FFT_PERIOD_MAX=280
GENE_FFT_RANKING_PERIODS="193,196,199"

RANDOMIZE_SEED=12345
RANDOMIZE_SEARCH_WINDOW=100000
RANDOMIZE_FALLBACK="uniform"

FRAG_COUNT_MIN=100
FRAG_COUNT_MAX=1000
FRAG_PLOT_MIN=100
FRAG_PLOT_MAX=1000
HEATMAP_MIN_FRAG=100
HEATMAP_MAX_FRAG=500
HEATMAP_NORMALIZATION="fragment-zscore"

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
  nucleosuite mnase-suite \
      --bam sample.bam \
      --fasta genome.fa \
      --ctcf-bed ctcf.bed \
      --outdir results \
      [options]

The same workflow can be run directly with mnase_full_suite.sh.

Inputs and execution:
  --bam FILE_OR_GLOB [MORE ...] One or more coordinate-sorted paired-end MNase BAMs.
                                Mutually exclusive with --fragments.
  --fragments FILE [MORE ...]   Fragment BED, BED.gz or bigBed files. Only the first three
                                columns are required. Mutually exclusive with --bam.
  --sample-name NAME            Output sample name. Default: derived from input filenames.
  --fasta FILE                  Matching indexed reference FASTA, or an indexable FASTA.
  --ctcf-bed FILE               BED3+ CTCF coordinates for aggregation and region extraction.
                                Optional with --resource-set hg19-gm12878.
  --outdir DIR                  Output directory.
  --cores N                     Process up to N contigs concurrently. Per-contig outputs
                                are written under OUTDIR/per_contig and combined under
                                OUTDIR/combined. Default: 1.
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

Regional and bundled resources:
  --states-bed FILE             BED3+ chromatin-state segmentation. Used for peak-distance
                                stratification, fragment-length profiles and gene sets.
  --genes-bed FILE              BED3+ gene regions, one unique record per gene. Requires
                                --states-bed and enables pooled gene-category DAC analyses.
  --gene-set-config FILE        TSV with set_name and include_rule columns and optional
                                exclude_if_candidate exclusions. If omitted, the bundled
                                active, weak and repressed classification is used.
  --expression FILE             Long-format expression TSV. When supplied, expression
                                analysis runs automatically in the combined workflow.
  --expression-value-column N   Expression value column name. Default: nTPM.
  --expression-gene-column N    Ensembl gene ID column name. Default: Gene.
  --expression-name-column N    Gene-name column name. Default: Gene name.
  --expression-profile-column N Expression profile/cell-line column. Default: Cell line.
  --expression-focus-profile N  Profile highlighted in plots; may be repeated.
  --tss-expression-resource F   Tissue-expression TSV/TSV.gz for TSS quintiles. Default: bundled HPA tissue consensus.
  --tss-expression-tissue N     Tissue/profile selector; use underscores for spaces. Default: bone_marrow.
  --tss-expression-window N     Bases on each side of TSS. Default: 2000.
  --venn-sets NAMES             Comma-separated two or three candidate sets for the Venn
                                diagram. Default: active_genes,weak_genes,repressed_genes.
  --states-label-column N       One-based state label column. Default: 4.
  --resource-set NAME           Bundled resource collection. Available: hg19-gm12878.
                                Supplies hg19 genes, GM12878 states, CTCF sites and default
                                gene-set rules unless the corresponding options are given.
  --interval-format VALUE       bed, bigbed or both for interval outputs. Default: both.
  --blacklist-bed FILE          Override the assembly-specific blacklist.
  --no-blacklist                Disable blacklist filtering. hg19 v2 is otherwise
                                enabled automatically for hg19/GRCh37 inputs.

Fragment selection:
  --pns-frag-lower N            PNS lower fragment length. Default: 120.
  --pns-frag-upper N            PNS upper fragment length. Default: 180.
  --pns-mode-length N           PNS modal fragment length. Default: 147.
  --bigbed-score-scale N        PNS peak multiplier used for integer bigBed scores. Default: 1000.
  --fine-frag-lower N           Ranged dyad/WW lower length. Default: 146.
  --fine-frag-upper N           Ranged dyad/WW upper length. Default: 148.
  --exact-size N                Exact dyad and fragment-end length. Default: 147.
  --max-duplicates N            Identical-fragment copy limit; 0 disables. Default: 1.
  --max-per-coordinate N        Optional dyad/end coordinate cap; 0 disables. Default: 0.
  --dedup-scope VALUE           all_bams or per_bam. Default: all_bams.
  --even-dyad VALUE             split, left or right. Default: split.

PNS and peak settings:
  --pns-smooth-window N         Optional Savitzky-Golay window; 0 disables. Default: 0.
  --pns-smooth-order N          Optional smoothing polynomial order. Default: 2.
  --pns-max-neg-run N           Zero-or-negative bases bridged within PNS peaks. Default: 0.
  --peak-smooth-window N        Standalone PNS peak-caller window; 0 disables. Default: 0.
  --peak-smooth-order N         Standalone peak-caller order. Default: 2.

DAC, NRL and distance settings:
  --dac-dmax N                  Maximum DAC distance. Default: 2000.
  --dac-window-size N           DAC genomic window size. Default: 100000.
  --dac-algorithm VALUE         auto, sparse or fft. Default: auto.
  --nrl-min-distance N          Main NRL lower bound. Default: 1.
  --nrl-max-distance N          Main NRL upper bound. Default: 1500.
  --nrl-peak-resolution N       Long-range NRL peak resolution in bp. Default: 160.
                                Detection smoothing uses resolution/3 and local-max smoothing uses resolution/6, snapped down to 10n+1 windows.
  --distance-x-major-tick N     Major x-axis tick interval for numeric distance plots. Default: automatic.
  --distance-x-minor-tick N     Minor x-axis tick interval. Default: derived from the major interval.
  --short-periodicity-min N     Short-periodicity lower bound. Default: 1.
  --short-periodicity-max N     Short-periodicity upper bound. Default: 144.
  --intermediate-periodicity-min N  Nucleosome-scale lower bound. Default: 150.
  --intermediate-periodicity-max N  Nucleosome-scale upper bound. Default: 220.
  --intermediate-periodicity-resolution N Peak resolution. Default: 8.
  --distance-adjacent-max N     Adjacent PNS peak-distance maximum. Default: 500.
  --distance-long-max N         1-7-order PNS distance maximum. Default: 1500.
  --distance-long-max-order N   Long-range neighbour orders. Default: 7.
  --state-distance-max N        ChromHMM overlay maximum adjacent distance. Default: 500.
  --state-distance-smooth-window N  State overlay Savitzky-Golay window. Default: 21.
  --state-distance-smooth-order N   State overlay polynomial order. Default: 2.
  --position-percentile-interval N  Width of independently ranked peak-score
                                    percentile groups used in both directional
                                    comparisons. Default: 25.
  --score-z-limit N              Symmetric score-correlation z-score axis limit;
                                    0 disables. Default: 10.
  --distance-histogram-x-max N   Displayed distance-histogram maximum in bp. Default: 300.
  --percentile-boxplot-y-max N   Displayed percentile-boxplot maximum in bp;
                                    0 disables. Default: 500.

Positive-run settings:
  --positive-runs-threshold N      A base is positive only when its score is greater than N. Default: 0.
  --positive-runs-chunk-size N     BigWig scan chunk size in bp. Default: 1000000.
  --positive-runs-min-length N     Minimum retained positive run length. Default: 1.
  --positive-runs-max-length N     Maximum retained run length; 0 disables. Default: 0.
  --positive-runs-plot-x-max N     Displayed run-length maximum; 0 uses all lengths. Default: 550.
  --positive-runs-normalization V  count, fraction or percent. Default: count.

Peak-score frequency settings:
  --peak-score-normalization V    count, fraction, percent or density. Default: count.

Gene-expression FFT settings (used only with --expression):
  --gene-fft-window N           Strand-aware window from each TSS. Default: 10000.
  --gene-fft-period-min N       Minimum FFT period. Default: 120.
  --gene-fft-period-max N       Maximum FFT period. Default: 280.
  --gene-fft-ranking-periods N  Comma-separated ranking periods. Default: 193,196,199.

Alignment and regional settings:
  --contigs VALUE [VALUE ...]   One or more contig selectors. Supports comma lists,
                                numeric ranges, autosomes, all, and genomic intervals.
                                Examples: chr20; 1,2,3; chr1-22 chrX chrY;
                                autosomes; all; chr20:1000000-2000000.
                                Default: autosomes.
  --ctcf-flank N                Region-extraction half-width. Default: 2000.
  --aggregate-window-half N     Aggregation half-width. Default: 2500.
  --region-peak-flank N         Peak search flank for region extraction. Default: 2000.

Fragment-length settings:
  --frag-count-min N            Count minimum. Default: 100.
  --frag-count-max N            Count maximum. Default: 1000.
  --frag-plot-min N             Plot minimum. Default: 100.
  --frag-plot-max N             Plot upper limit; stops at longest count. Default: 1000.
  --heatmap-min-frag N          Heatmap minimum. Default: 100.
  --heatmap-max-frag N          Heatmap maximum. Default: 500.
  --heatmap-normalization NAME  fragment-zscore, profile-percent,
                                fragment-percent, profile-minmax,
                                fragment-minmax or none. Default: fragment-zscore.

Randomized-control mode:
  --randomize                   Run a randomized control instead of the observed analysis.
                                Randomization occurs in 00_setup and the complete normal
                                analysis tree is then generated with _randomized_control names.
  --randomize-seed N            Reproducible dinucleotide-matched randomization seed. Default: 12345.
  --randomize-search-window N   Local dinucleotide search window. Default: 100000.
  --randomize-fallback VALUE    uniform or skip. Default: uniform.

Optional skips:
  --skip-nrl
  --skip-region-extract
  --skip-fragment-heatmap
  --skip-tss-expression-quintiles Skip tissue-expression quintile TSS analysis.
  --skip-gene-expression
  --skip-positive-runs
  --skip-peak-score-frequency

Other:
  --force                       Re-run completed steps.
  --resume                      Reuse only outputs with matching parameters and completed markers.
  --dry-run                     Validate inputs and print planned commands without running them.
  -h, --help                    Show this help.

Environment overrides:
  NUCLEOSUITE_BIN               NucleoSuite executable. Default: nucleosuite.
  PYTHON_BIN                    Python executable. Default: python.
EOF
}

require_value() {
    [[ $# -ge 2 && -n "${2:-}" ]] || { echo "ERROR: $1 requires a value" >&2; exit 2; }
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)
            shift
            [[ $# -gt 0 && "$1" != -* ]] || { echo "ERROR: --bam requires at least one BAM path or pattern" >&2; exit 2; }
            while [[ $# -gt 0 && "$1" != -* ]]; do
                BAM_INPUTS+=("$1")
                shift
            done
            ;;
        --bam=*)
            [[ -n "${1#*=}" ]] || { echo "ERROR: --bam requires a BAM path or pattern" >&2; exit 2; }
            BAM_INPUTS+=("${1#*=}")
            shift
            ;;
        --fragments)
            shift
            [[ $# -gt 0 && "$1" != -* ]] || { echo "ERROR: --fragments requires at least one path or pattern" >&2; exit 2; }
            while [[ $# -gt 0 && "$1" != -* ]]; do
                FRAGMENT_INPUTS+=("$1")
                shift
            done
            ;;
        --fragments=*)
            [[ -n "${1#*=}" ]] || { echo "ERROR: --fragments requires a path or pattern" >&2; exit 2; }
            FRAGMENT_INPUTS+=("${1#*=}")
            shift
            ;;
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
        --states-label-column) require_value "$@"; STATES_LABEL_COLUMN="$2"; shift 2 ;;
        --resource-set) require_value "$@"; RESOURCE_SET="$2"; shift 2 ;;
        --outdir) require_value "$@"; OUTDIR="$2"; shift 2 ;;
        --combine-bigwig-method)
            require_value "$@"
            [[ "$2" == "direct" || "$2" == "bedgraph" || "$2" == "bedgraphs" ]] || { echo "ERROR: --combine-bigwig-method must be direct or bedgraph" >&2; exit 2; }
            shift 2
            ;; # consumed by the Python wrapper
        --combine-bigwig-method=*)
            WRAPPER_BIGWIG_METHOD="${1#*=}"
            [[ "$WRAPPER_BIGWIG_METHOD" == "direct" || "$WRAPPER_BIGWIG_METHOD" == "bedgraph" || "$WRAPPER_BIGWIG_METHOD" == "bedgraphs" ]] || { echo "ERROR: --combine-bigwig-method must be direct or bedgraph" >&2; exit 2; }
            shift
            ;; # consumed by the Python wrapper
        --analysis-scope)
            require_value "$@"
            WRAPPER_ANALYSIS_SCOPE="$2"
            [[ "$WRAPPER_ANALYSIS_SCOPE" == "combined-only" || "$WRAPPER_ANALYSIS_SCOPE" == "per-contig-and-combined" ]] || { echo "ERROR: invalid --analysis-scope" >&2; exit 2; }
            shift 2
            ;; # consumed by the Python wrapper
        --analysis-scope=*)
            WRAPPER_ANALYSIS_SCOPE="${1#*=}"
            [[ "$WRAPPER_ANALYSIS_SCOPE" == "combined-only" || "$WRAPPER_ANALYSIS_SCOPE" == "per-contig-and-combined" ]] || { echo "ERROR: invalid --analysis-scope" >&2; exit 2; }
            shift
            ;; # consumed by the Python wrapper
        --combine-prerequisites-only) COMBINE_PREREQUISITES_ONLY=1; shift ;;
        --cores)
            require_value "$@"
            [[ "$2" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --cores must be a positive integer" >&2; exit 2; }
            shift 2
            ;; # applied by the Python multicontig launcher
        --cores=*)
            WRAPPER_CORES="${1#*=}"
            [[ "$WRAPPER_CORES" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --cores must be a positive integer" >&2; exit 2; }
            shift
            ;; # applied by the Python multicontig launcher
        --contigs)
            shift
            [[ $# -gt 0 && "$1" != -* ]] || { echo "ERROR: --contigs requires at least one selector" >&2; exit 2; }
            CONTIGS=()
            while [[ $# -gt 0 && "$1" != -* ]]; do
                CONTIGS+=("$1")
                shift
            done
            ;;
        --contigs=*)
            [[ -n "${1#*=}" ]] || { echo "ERROR: --contigs requires a selector" >&2; exit 2; }
            CONTIGS=("${1#*=}")
            shift
            ;;
        --interval-format) require_value "$@"; INTERVAL_FORMAT="$2"; shift 2 ;;
        --pns-frag-lower) require_value "$@"; PNS_FRAG_LOWER="$2"; shift 2 ;;
        --pns-frag-upper) require_value "$@"; PNS_FRAG_UPPER="$2"; shift 2 ;;
        --broad-frag-lower) require_value "$@"; PNS_FRAG_LOWER="$2"; shift 2 ;;
        --broad-frag-upper) require_value "$@"; PNS_FRAG_UPPER="$2"; shift 2 ;;
        --pns-mode-length) require_value "$@"; PNS_MODE_LENGTH="$2"; shift 2 ;;
    --bigbed-score-scale) require_value "$@"; BIGBED_SCORE_SCALE="$2"; shift 2 ;;
        --fine-frag-lower) require_value "$@"; FINE_FRAG_LOWER="$2"; shift 2 ;;
        --fine-frag-upper) require_value "$@"; FINE_FRAG_UPPER="$2"; shift 2 ;;
        --exact-size) require_value "$@"; EXACT_SIZE="$2"; shift 2 ;;
        --max-duplicates) require_value "$@"; MAX_DUPLICATES="$2"; shift 2 ;;
        --max-per-coordinate) require_value "$@"; MAX_PER_COORDINATE="$2"; shift 2 ;;
        --dedup-scope) require_value "$@"; DEDUP_SCOPE="$2"; shift 2 ;;
        --even-dyad) require_value "$@"; EVEN_DYAD="$2"; shift 2 ;;
        --pns-smooth-window) require_value "$@"; PNS_SMOOTH_WINDOW="$2"; shift 2 ;;
        --pns-smooth-order) require_value "$@"; PNS_SMOOTH_ORDER="$2"; shift 2 ;;
        --pns-max-neg-run) require_value "$@"; PNS_MAX_NEG_RUN="$2"; shift 2 ;;
        --peak-smooth-window) require_value "$@"; PEAK_SMOOTH_WINDOW="$2"; shift 2 ;;
        --peak-smooth-order) require_value "$@"; PEAK_SMOOTH_ORDER="$2"; shift 2 ;;
        --ctcf-flank) require_value "$@"; CTCF_FLANK="$2"; shift 2 ;;
        --aggregate-window-half) require_value "$@"; AGGREGATE_WINDOW_HALF="$2"; shift 2 ;;
        --region-peak-flank) require_value "$@"; REGION_PEAK_FLANK="$2"; shift 2 ;;
        --dac-dmax) require_value "$@"; DAC_DMAX="$2"; shift 2 ;;
        --dac-window-size) require_value "$@"; DAC_WINDOW_SIZE="$2"; shift 2 ;;
        --dac-algorithm) require_value "$@"; DAC_ALGORITHM="$2"; shift 2 ;;
        --nrl-min-distance) require_value "$@"; NRL_MIN_DISTANCE="$2"; shift 2 ;;
        --nrl-max-distance) require_value "$@"; NRL_MAX_DISTANCE="$2"; shift 2 ;;
        --nrl-peak-resolution) require_value "$@"; NRL_PEAK_RESOLUTION="$2"; shift 2 ;;
        --distance-x-major-tick) require_value "$@"; DISTANCE_X_MAJOR_TICK="$2"; shift 2 ;;
        --distance-x-minor-tick) require_value "$@"; DISTANCE_X_MINOR_TICK="$2"; shift 2 ;;
        --short-periodicity-min) require_value "$@"; SHORT_PERIODICITY_MIN="$2"; shift 2 ;;
        --short-periodicity-max) require_value "$@"; SHORT_PERIODICITY_MAX="$2"; shift 2 ;;
        --intermediate-periodicity-min) require_value "$@"; INTERMEDIATE_PERIODICITY_MIN="$2"; shift 2 ;;
        --intermediate-periodicity-max) require_value "$@"; INTERMEDIATE_PERIODICITY_MAX="$2"; shift 2 ;;
        --intermediate-periodicity-resolution) require_value "$@"; INTERMEDIATE_PERIODICITY_RESOLUTION="$2"; shift 2 ;;
        --distance-adjacent-max) require_value "$@"; DISTANCE_ADJACENT_MAX="$2"; shift 2 ;;
        --distance-long-max) require_value "$@"; DISTANCE_LONG_MAX="$2"; shift 2 ;;
        --distance-long-max-order) require_value "$@"; DISTANCE_LONG_MAX_ORDER="$2"; shift 2 ;;
        --state-distance-max) require_value "$@"; STATE_DISTANCE_MAX="$2"; shift 2 ;;
        --state-distance-smooth-window) require_value "$@"; STATE_DISTANCE_SMOOTH_WINDOW="$2"; shift 2 ;;
        --state-distance-smooth-order) require_value "$@"; STATE_DISTANCE_SMOOTH_ORDER="$2"; shift 2 ;;
        --position-percentile-interval) require_value "$@"; POSITION_PERCENTILE_INTERVAL="$2"; shift 2 ;;
        --score-z-limit) require_value "$@"; SCORE_Z_LIMIT="$2"; shift 2 ;;
        --distance-histogram-x-max) require_value "$@"; DISTANCE_HISTOGRAM_X_MAX="$2"; shift 2 ;;
        --percentile-boxplot-y-max) require_value "$@"; PERCENTILE_BOXPLOT_Y_MAX="$2"; shift 2 ;;
        --positive-runs-threshold) require_value "$@"; POSITIVE_RUNS_THRESHOLD="$2"; shift 2 ;;
        --positive-runs-chunk-size) require_value "$@"; POSITIVE_RUNS_CHUNK_SIZE="$2"; shift 2 ;;
        --positive-runs-min-length) require_value "$@"; POSITIVE_RUNS_MIN_LENGTH="$2"; shift 2 ;;
        --positive-runs-max-length) require_value "$@"; POSITIVE_RUNS_MAX_LENGTH="$2"; shift 2 ;;
        --positive-runs-plot-x-max) require_value "$@"; POSITIVE_RUNS_PLOT_X_MAX="$2"; shift 2 ;;
        --positive-runs-normalization) require_value "$@"; POSITIVE_RUNS_NORMALIZATION="$2"; shift 2 ;;
        --peak-score-normalization) require_value "$@"; PEAK_SCORE_NORMALIZATION="$2"; shift 2 ;;
        --gene-fft-window) require_value "$@"; GENE_FFT_WINDOW="$2"; shift 2 ;;
        --gene-fft-period-min) require_value "$@"; GENE_FFT_PERIOD_MIN="$2"; shift 2 ;;
        --gene-fft-period-max) require_value "$@"; GENE_FFT_PERIOD_MAX="$2"; shift 2 ;;
        --gene-fft-ranking-periods) require_value "$@"; GENE_FFT_RANKING_PERIODS="$2"; shift 2 ;;
        --randomize-seed) require_value "$@"; RANDOMIZE_SEED="$2"; shift 2 ;;
        --randomize-search-window) require_value "$@"; RANDOMIZE_SEARCH_WINDOW="$2"; shift 2 ;;
        --randomize-fallback) require_value "$@"; RANDOMIZE_FALLBACK="$2"; shift 2 ;;
        --frag-count-min) require_value "$@"; FRAG_COUNT_MIN="$2"; shift 2 ;;
        --frag-count-max) require_value "$@"; FRAG_COUNT_MAX="$2"; shift 2 ;;
        --frag-plot-min) require_value "$@"; FRAG_PLOT_MIN="$2"; shift 2 ;;
        --frag-plot-max) require_value "$@"; FRAG_PLOT_MAX="$2"; shift 2 ;;
        --heatmap-min-frag) require_value "$@"; HEATMAP_MIN_FRAG="$2"; shift 2 ;;
        --heatmap-max-frag) require_value "$@"; HEATMAP_MAX_FRAG="$2"; shift 2 ;;
        --heatmap-normalization) require_value "$@"; HEATMAP_NORMALIZATION="$2"; shift 2 ;;
        --skip-nrl) SKIP_NRL=1; shift ;;
        --skip-region-extract) SKIP_REGION_EXTRACT=1; shift ;;
        --skip-fragment-heatmap) SKIP_FRAGMENT_HEATMAP=1; shift ;;
        --skip-gene-expression) SKIP_GENE_EXPRESSION=1; shift ;;
        --skip-tss-expression-quintiles) SKIP_TSS_EXPRESSION_QUINTILES=1; shift ;;
        --randomize) RUN_MODE="randomized"; shift ;;
        --randomized-control-input) RUN_MODE="randomized"; RANDOMIZED_INPUT_READY=1; shift ;;
        --provenance-bam) require_value "$@"; PROVENANCE_BAMS+=("$2"); shift 2 ;;
        --provenance-fragment) require_value "$@"; PROVENANCE_FRAGMENTS+=("$2"); shift 2 ;;
        --trusted-combined-prerequisites) TRUST_EXISTING_OUTPUTS=1; shift ;;
        --skip-positive-runs) SKIP_POSITIVE_RUNS=1; shift ;;
        --skip-peak-score-frequency) SKIP_PEAK_SCORE_FREQUENCY=1; shift ;;
        --force) FORCE=1; shift ;;
        --resume) REUSE_EXISTING_OUTPUTS=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --validate-only) VALIDATE_ONLY=1; shift ;;
        -h|--help) usage; exit 0 ;;
    --help-plotting) usage; plotting_usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

fatal() { echo "ERROR: $*" >&2; exit 1; }
check_uint() { [[ "$2" =~ ^[0-9]+$ ]] || fatal "$1 must be a non-negative integer"; }
check_posint() { [[ "$2" =~ ^[1-9][0-9]*$ ]] || fatal "$1 must be a positive integer"; }
check_nonneg_number() { [[ "$2" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] || fatal "$1 must be a non-negative number"; }
check_pos_number() { [[ "$2" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && awk -v v="$2" 'BEGIN{exit !(v>0)}' || fatal "$1 must be a positive number"; }

if [[ ${#BAM_INPUTS[@]} -gt 0 && ${#FRAGMENT_INPUTS[@]} -gt 0 ]]; then fatal "Use either --bam or --fragments, not both"; fi
if [[ ${#BAM_INPUTS[@]} -eq 0 && ${#FRAGMENT_INPUTS[@]} -eq 0 ]]; then fatal "One of --bam or --fragments is required"; fi
[[ -n "$FASTA" ]] || fatal "--fasta is required"
[[ -n "$OUTDIR" ]] || fatal "--outdir is required"
[[ -f "$FASTA" ]] || fatal "FASTA not found: $FASTA"
[[ -z "$BLACKLIST_BED" || -f "$BLACKLIST_BED" ]] || fatal "Blacklist BED not found: $BLACKLIST_BED"

expand_bam_inputs() {
    local token match
    local -a matches=() expanded=()
    for token in "${BAM_INPUTS[@]}"; do
        matches=()
        if [[ -f "$token" ]]; then
            matches=("$token")
        else
            mapfile -t matches < <(compgen -G "$token" || true)
        fi
        [[ ${#matches[@]} -gt 0 ]] || fatal "BAM input did not match any files: $token"
        for match in "${matches[@]}"; do
            [[ -f "$match" ]] || continue
            [[ "${match,,}" == *.bam ]] || fatal "BAM input is not a .bam file: $match"
            expanded+=("$(realpath "$match")")
        done
    done
    [[ ${#expanded[@]} -gt 0 ]] || fatal "No BAM files were found"
    mapfile -t BAMS < <(printf '%s\n' "${expanded[@]}" | sort -V -u)
}
if [[ ${#BAM_INPUTS[@]} -gt 0 ]]; then
    expand_bam_inputs
    INPUT_MODE="bam"
else
    for token in "${FRAGMENT_INPUTS[@]}"; do
        if [[ -f "$token" ]]; then
            FRAGMENTS+=("$(realpath "$token")")
        else
            mapfile -t matches < <(compgen -G "$token" || true)
            [[ ${#matches[@]} -gt 0 ]] || fatal "Fragment input did not match any files: $token"
            for match in "${matches[@]}"; do FRAGMENTS+=("$(realpath "$match")"); done
        fi
    done
    mapfile -t FRAGMENTS < <(printf '%s\n' "${FRAGMENTS[@]}" | sort -V -u)
    INPUT_MODE="fragments"
fi
if [[ "${#PROVENANCE_BAMS[@]}" -eq 0 && "${#PROVENANCE_FRAGMENTS[@]}" -eq 0 ]]; then
    if [[ "$INPUT_MODE" == "bam" ]]; then
        PROVENANCE_BAMS=("${BAMS[@]}")
    else
        PROVENANCE_FRAGMENTS=("${FRAGMENTS[@]}")
    fi
fi

BLACKLIST_ARGS=()
[[ -n "$BLACKLIST_BED" ]] && BLACKLIST_ARGS=(--blacklist-bed "$BLACKLIST_BED")
ANALYSIS_INPUT_ARGS=()
SOURCE_REFERENCE_ARGS=()
if [[ "$INPUT_MODE" == "bam" ]]; then
    ANALYSIS_INPUT_ARGS=(-b "${BAMS[@]}")
else
    ANALYSIS_INPUT_ARGS=(--fragments "${FRAGMENTS[@]}")
    SOURCE_REFERENCE_ARGS=(--fasta "$FASTA")
fi
SOURCE_INPUT_MODE="$INPUT_MODE"
SOURCE_INPUT_ARGS=("${ANALYSIS_INPUT_ARGS[@]}")

derive_sample_name() {
    if [[ -n "$SAMPLE_NAME" ]]; then
        SAMPLE_NAME="${SAMPLE_NAME//[^A-Za-z0-9._-]/_}"
        [[ -n "$SAMPLE_NAME" ]] || fatal "--sample-name does not contain any usable characters"
        return
    fi
    local path stem stripped candidate=""
    local consistent=1
    shopt -s nocasematch
    for path in "${BAMS[@]}"; do
        stem="$(basename "$path")"
        stem="${stem%.bam}"
        stripped="$stem"
        if [[ "$stem" =~ ^(.*)[._-]chr([0-9]+|X|Y|M|MT)$ ]]; then
            stripped="${BASH_REMATCH[1]}"
        fi
        if [[ -z "$candidate" ]]; then
            candidate="$stripped"
        elif [[ "$candidate" != "$stripped" ]]; then
            consistent=0
        fi
    done
    shopt -u nocasematch
    if [[ "$INPUT_MODE" == "bam" && ${#BAMS[@]} -eq 1 ]]; then
        candidate="$(basename "${BAMS[0]}")"
        candidate="${candidate%.bam}"
    elif [[ "$INPUT_MODE" == "fragments" && ${#FRAGMENTS[@]} -eq 1 ]]; then
        candidate="$(basename "${FRAGMENTS[0]}")"
        candidate="${candidate%.bed.gz}"
        candidate="${candidate%.bed}"
        candidate="${candidate%.bb}"
        candidate="${candidate%.bigBed}"
    elif [[ "$consistent" -eq 0 || -z "$candidate" ]]; then
        candidate="multi_${INPUT_MODE}"
    fi
    SAMPLE_NAME="${candidate//[^A-Za-z0-9._-]/_}"
}
derive_sample_name

for pair in \
    "--pns-frag-lower:$PNS_FRAG_LOWER" "--pns-frag-upper:$PNS_FRAG_UPPER" \
    "--pns-mode-length:$PNS_MODE_LENGTH" "--fine-frag-lower:$FINE_FRAG_LOWER" \
    "--fine-frag-upper:$FINE_FRAG_UPPER" "--exact-size:$EXACT_SIZE" \
    "--dac-dmax:$DAC_DMAX" \
    "--nrl-min-distance:$NRL_MIN_DISTANCE" "--nrl-max-distance:$NRL_MAX_DISTANCE" \
    "--nrl-peak-resolution:$NRL_PEAK_RESOLUTION" \
    "--short-periodicity-min:$SHORT_PERIODICITY_MIN" "--short-periodicity-max:$SHORT_PERIODICITY_MAX" \
    "--intermediate-periodicity-min:$INTERMEDIATE_PERIODICITY_MIN" "--intermediate-periodicity-max:$INTERMEDIATE_PERIODICITY_MAX" "--intermediate-periodicity-resolution:$INTERMEDIATE_PERIODICITY_RESOLUTION" \
    "--aggregate-window-half:$AGGREGATE_WINDOW_HALF" "--state-distance-max:$STATE_DISTANCE_MAX" \
    "--state-distance-smooth-window:$STATE_DISTANCE_SMOOTH_WINDOW" "--position-percentile-interval:$POSITION_PERCENTILE_INTERVAL" \
    "--gene-fft-window:$GENE_FFT_WINDOW" \
    "--gene-fft-period-min:$GENE_FFT_PERIOD_MIN" "--gene-fft-period-max:$GENE_FFT_PERIOD_MAX" \
    "--randomize-search-window:$RANDOMIZE_SEARCH_WINDOW"; do
    check_posint "${pair%%:*}" "${pair#*:}"
done
for pair in \
    "--max-duplicates:$MAX_DUPLICATES" "--max-per-coordinate:$MAX_PER_COORDINATE" "--ctcf-flank:$CTCF_FLANK" \
    "--region-peak-flank:$REGION_PEAK_FLANK"; do
    check_uint "${pair%%:*}" "${pair#*:}"
done
for pair in \
    "--pns-smooth-window:$PNS_SMOOTH_WINDOW" "--pns-smooth-order:$PNS_SMOOTH_ORDER" \
    "--pns-max-neg-run:$PNS_MAX_NEG_RUN" "--peak-smooth-window:$PEAK_SMOOTH_WINDOW" \
    "--peak-smooth-order:$PEAK_SMOOTH_ORDER"; do
    check_uint "${pair%%:*}" "${pair#*:}"
done
if [[ "$PNS_SMOOTH_WINDOW" -gt 0 ]]; then
    (( PNS_SMOOTH_WINDOW >= 3 && PNS_SMOOTH_WINDOW % 2 == 1 )) || fatal "--pns-smooth-window must be 0 or an odd integer of at least 3"
    (( PNS_SMOOTH_ORDER < PNS_SMOOTH_WINDOW )) || fatal "--pns-smooth-order must be smaller than --pns-smooth-window"
fi
if [[ "$PEAK_SMOOTH_WINDOW" -gt 0 ]]; then
    (( PEAK_SMOOTH_WINDOW >= 3 && PEAK_SMOOTH_WINDOW % 2 == 1 )) || fatal "--peak-smooth-window must be 0 or an odd integer of at least 3"
    (( PEAK_SMOOTH_ORDER < PEAK_SMOOTH_WINDOW )) || fatal "--peak-smooth-order must be smaller than --peak-smooth-window"
fi
(( PNS_FRAG_LOWER <= PNS_FRAG_UPPER )) || fatal "PNS fragment lower exceeds upper"
(( FINE_FRAG_LOWER <= FINE_FRAG_UPPER )) || fatal "fine fragment lower exceeds upper"
(( FRAG_COUNT_MIN <= FRAG_COUNT_MAX )) || fatal "fragment count minimum exceeds maximum"
(( FRAG_PLOT_MIN <= FRAG_PLOT_MAX )) || fatal "fragment plot minimum exceeds maximum"
(( HEATMAP_MIN_FRAG <= HEATMAP_MAX_FRAG )) || fatal "heatmap minimum exceeds maximum"
(( NRL_MIN_DISTANCE < NRL_MAX_DISTANCE )) || fatal "NRL minimum distance must be less than maximum"
[[ -z "$DISTANCE_X_MAJOR_TICK" ]] || check_pos_number "--distance-x-major-tick" "$DISTANCE_X_MAJOR_TICK"
[[ -z "$DISTANCE_X_MINOR_TICK" ]] || check_pos_number "--distance-x-minor-tick" "$DISTANCE_X_MINOR_TICK"
(( SHORT_PERIODICITY_MIN < SHORT_PERIODICITY_MAX )) || fatal "short-periodicity minimum must be less than maximum"
(( INTERMEDIATE_PERIODICITY_MIN < INTERMEDIATE_PERIODICITY_MAX )) || fatal "intermediate-periodicity minimum must be less than maximum"
(( GENE_FFT_PERIOD_MIN < GENE_FFT_PERIOD_MAX )) || fatal "gene FFT minimum period must be less than maximum"
(( STATE_DISTANCE_SMOOTH_WINDOW % 2 == 1 )) || fatal "--state-distance-smooth-window must be odd"
check_uint "--state-distance-smooth-order" "$STATE_DISTANCE_SMOOTH_ORDER"
(( STATE_DISTANCE_SMOOTH_ORDER < STATE_DISTANCE_SMOOTH_WINDOW )) || fatal "state-distance smoothing order must be smaller than the window"
(( POSITION_PERCENTILE_INTERVAL <= 100 )) || fatal "--position-percentile-interval must be between 1 and 100"
check_nonneg_number "--score-z-limit" "$SCORE_Z_LIMIT"
check_pos_number "--distance-histogram-x-max" "$DISTANCE_HISTOGRAM_X_MAX"
check_nonneg_number "--percentile-boxplot-y-max" "$PERCENTILE_BOXPLOT_Y_MAX"
check_nonneg_number "--positive-runs-threshold" "$POSITIVE_RUNS_THRESHOLD"
check_posint "--positive-runs-chunk-size" "$POSITIVE_RUNS_CHUNK_SIZE"
check_posint "--positive-runs-min-length" "$POSITIVE_RUNS_MIN_LENGTH"
check_uint "--positive-runs-max-length" "$POSITIVE_RUNS_MAX_LENGTH"
check_uint "--positive-runs-plot-x-max" "$POSITIVE_RUNS_PLOT_X_MAX"
if [[ "$POSITIVE_RUNS_MAX_LENGTH" -gt 0 && "$POSITIVE_RUNS_MAX_LENGTH" -lt "$POSITIVE_RUNS_MIN_LENGTH" ]]; then fatal "--positive-runs-max-length cannot be smaller than --positive-runs-min-length"; fi
[[ "$POSITIVE_RUNS_NORMALIZATION" =~ ^(count|fraction|percent)$ ]] || fatal "--positive-runs-normalization must be count, fraction or percent"
[[ "$PEAK_SCORE_NORMALIZATION" =~ ^(count|fraction|percent|density)$ ]] || fatal "--peak-score-normalization must be count, fraction, percent or density"
if [[ "$SKIP_NRL" -eq 0 ]]; then
    required_dac_distance=$NRL_MAX_DISTANCE
    (( SHORT_PERIODICITY_MAX > required_dac_distance )) && required_dac_distance=$SHORT_PERIODICITY_MAX
    (( INTERMEDIATE_PERIODICITY_MAX > required_dac_distance )) && required_dac_distance=$INTERMEDIATE_PERIODICITY_MAX
    (( DAC_DMAX >= required_dac_distance )) || fatal "--dac-dmax must be at least ${required_dac_distance} for the requested NRL analyses"
fi
[[ "$MAX_DUPLICATES" =~ ^[0-9]+$ ]] || fatal "--max-duplicates must be a non-negative integer"
[[ "$MAX_PER_COORDINATE" =~ ^[0-9]+$ ]] || fatal "--max-per-coordinate must be a non-negative integer"
[[ "$DEDUP_SCOPE" == "all_bams" || "$DEDUP_SCOPE" == "per_bam" ]] || fatal "invalid --dedup-scope"
[[ "$EVEN_DYAD" == "split" || "$EVEN_DYAD" == "left" || "$EVEN_DYAD" == "right" ]] || fatal "invalid --even-dyad"
[[ "$DAC_ALGORITHM" == "auto" || "$DAC_ALGORITHM" == "sparse" || "$DAC_ALGORITHM" == "fft" ]] || fatal "invalid --dac-algorithm"
[[ "$RANDOMIZE_FALLBACK" == "uniform" || "$RANDOMIZE_FALLBACK" == "skip" ]] || fatal "invalid --randomize-fallback"
check_uint "--randomize-seed" "$RANDOMIZE_SEED"
check_posint "--states-label-column" "$STATES_LABEL_COLUMN"

command -v "$NUCLEOSUITE_BIN" >/dev/null 2>&1 || fatal "NucleoSuite executable not found: $NUCLEOSUITE_BIN"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || fatal "Python executable not found: $PYTHON_BIN"
[[ "$INTERVAL_FORMAT" =~ ^(bed|bigbed|both)$ ]] || fatal "--interval-format must be bed, bigbed or both"
if [[ "$INTERVAL_FORMAT" != "bed" ]]; then
    command -v bedToBigBed >/dev/null 2>&1 || fatal "bedToBigBed is required for $INTERVAL_FORMAT interval output"
    command -v bigBedToBed >/dev/null 2>&1 || fatal "bigBedToBed is required for downstream bigBed input"
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
if [[ -n "$GENE_SET_CONFIG" && ( -z "$GENES_BED" || -z "$STATES_BED" ) ]]; then fatal "--gene-set-config requires both --genes-bed and --states-bed"; fi
if [[ -n "$STATES_BED" && -n "$GENES_BED" && -z "$GENE_SET_CONFIG" ]]; then
    GENE_SET_CONFIG="$($NUCLEOSUITE_BIN resources path default-gene-sets)"
fi
if [[ -n "$EXPRESSION" && -z "$GENES_BED" ]]; then fatal "--expression requires --genes-bed or a resource set that supplies genes"; fi

FASTA="$(realpath "$FASTA")"
CTCF_BED="$(realpath "$CTCF_BED")"
[[ -z "$STATES_BED" ]] || STATES_BED="$(realpath "$STATES_BED")"
[[ -z "$GENES_BED" ]] || GENES_BED="$(realpath "$GENES_BED")"
[[ -z "$GENE_SET_CONFIG" ]] || GENE_SET_CONFIG="$(realpath "$GENE_SET_CONFIG")"
[[ -z "$EXPRESSION" ]] || EXPRESSION="$(realpath "$EXPRESSION")"
TSS_EXPRESSION_RESOURCE="$(realpath "$TSS_EXPRESSION_RESOURCE")"
if [[ "$VALIDATE_ONLY" -eq 1 || "$DRY_RUN" -eq 1 ]]; then
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf 'NucleoSuite MNase suite dry run\nmode\t%s\noutdir\t%s\nblacklist\t%s\n' \
      "$RUN_MODE" "$OUTDIR" "${BLACKLIST_BED:-auto-if-hg19}"
    printf 'stages\tsetup,tracks,dac,nrl,scaling,aggregates,distances,region-extract,fragment-lengths,heatmaps,gene-expression,positive-runs,peak-analysis\n'
  fi
  exit 0
fi

mkdir -p "$OUTDIR"
OUTDIR="$(realpath "$OUTDIR")"

SAMPLE="$SAMPLE_NAME"
if [[ "$RUN_MODE" == "randomized" && "$SAMPLE" != *_randomized_control* ]]; then
    SAMPLE="${SAMPLE}_randomized_control"
    SAMPLE_NAME="$SAMPLE"
fi
SUPPORT_PREFIX=""
if [[ "$RUN_MODE" == "randomized" ]]; then
    SUPPORT_PREFIX="${SAMPLE}_"
fi
SETUP_DIR="$OUTDIR/00_setup"
GENE_SET_DIR="$OUTDIR/00_gene_sets"
COMBINED_TRACK_DIR="$OUTDIR/01_combined_tracks"
PNS_DIR="$COMBINED_TRACK_DIR/pns"
DYAD_EXACT_DIR="$COMBINED_TRACK_DIR/dyads/exact"
DYAD_RANGE_DIR="$COMBINED_TRACK_DIR/dyads/ranges"
ENDS_EXACT_DIR="$COMBINED_TRACK_DIR/fragment_ends/exact"
ENDS_RANGE_DIR="$COMBINED_TRACK_DIR/fragment_ends/ranges"
SEQUENCE_DIR="$COMBINED_TRACK_DIR/sequence"
DAC_DIR="$OUTDIR/02_dac"
NRL_DIR="$OUTDIR/04_nrl"
CTCF_AGG_DIR="$OUTDIR/05_ctcf_aggregation"
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
SCALED_DIR="$COMBINED_TRACK_DIR/scaled"
LOG_DIR="$OUTDIR/logs"
DONE_DIR="$OUTDIR/.done"
mkdir -p "$SETUP_DIR" "$GENE_SET_DIR" "$COMBINED_TRACK_DIR" "$PNS_DIR" "$SCALED_DIR" \
    "$DYAD_EXACT_DIR" "$DYAD_RANGE_DIR" "$ENDS_EXACT_DIR" "$ENDS_RANGE_DIR" "$SEQUENCE_DIR" \
    "$DAC_DIR" "$NRL_DIR" "$CTCF_AGG_DIR" "$TSS_AGG_DIR" "$DIST_DIR" "$REGION_DIR" \
    "$FRAG_DIR" "$HEATMAP_DIR" "$GENE_EXPRESSION_DIR" "$POSITIVE_RUNS_DIR" "$PEAK_ANALYSIS_DIR" \
    "$PEAK_SCORE_DIR" "$LOG_DIR" "$DONE_DIR"

materialize_bigbed_input() {
    local source="$1" label="$2"
    if [[ "${source,,}" == *.bb || "${source,,}" == *.bigbed ]]; then
        local target="$SETUP_DIR/${SUPPORT_PREFIX}${label}.bed"
        bigBedToBed "$source" "$target"
        printf '%s\n' "$target"
    else
        printf '%s\n' "$source"
    fi
}
CTCF_BED="$(materialize_bigbed_input "$CTCF_BED" ctcf_input)"
[[ -z "$STATES_BED" ]] || STATES_BED="$(materialize_bigbed_input "$STATES_BED" states_input)"
[[ -z "$GENES_BED" ]] || GENES_BED="$(materialize_bigbed_input "$GENES_BED" genes_input)"
INTERVAL_EXT="bed"
[[ "$INTERVAL_FORMAT" == "bigbed" ]] && INTERVAL_EXT="bb"

PARAMETERS="$SETUP_DIR/${SUPPORT_PREFIX}run_parameters.tsv"
{
    echo -e "parameter\tvalue"
    for name in SAMPLE_NAME RUN_MODE FASTA BLACKLIST_BED NO_BLACKLIST CTCF_BED STATES_BED GENES_BED GENE_SET_CONFIG EXPRESSION EXPRESSION_VALUE_COLUMN EXPRESSION_GENE_COLUMN EXPRESSION_NAME_COLUMN EXPRESSION_PROFILE_COLUMN TSS_EXPRESSION_RESOURCE TSS_EXPRESSION_TISSUE TSS_EXPRESSION_WINDOW RESOURCE_SET VENN_SETS OUTDIR INTERVAL_FORMAT \
        PNS_FRAG_LOWER PNS_FRAG_UPPER PNS_MODE_LENGTH FINE_FRAG_LOWER FINE_FRAG_UPPER \
        EXACT_SIZE DINUC_EXACT_A DINUC_EXACT_B MAX_DUPLICATES MAX_PER_COORDINATE DEDUP_SCOPE EVEN_DYAD PNS_SMOOTH_WINDOW \
        PNS_SMOOTH_ORDER PNS_MAX_NEG_RUN \
        PEAK_SMOOTH_WINDOW PEAK_SMOOTH_ORDER CTCF_FLANK AGGREGATE_WINDOW_HALF \
        REGION_PEAK_FLANK STATES_LABEL_COLUMN DAC_DMAX DAC_WINDOW_SIZE \
        DAC_ALGORITHM NRL_MIN_DISTANCE NRL_MAX_DISTANCE NRL_PEAK_RESOLUTION \
        SHORT_PERIODICITY_MIN SHORT_PERIODICITY_MAX INTERMEDIATE_PERIODICITY_MIN INTERMEDIATE_PERIODICITY_MAX INTERMEDIATE_PERIODICITY_RESOLUTION \
        DISTANCE_ADJACENT_MAX DISTANCE_LONG_MAX DISTANCE_LONG_MAX_ORDER STATE_DISTANCE_MAX STATE_DISTANCE_SMOOTH_WINDOW STATE_DISTANCE_SMOOTH_ORDER POSITION_PERCENTILE_INTERVAL SCORE_Z_LIMIT DISTANCE_HISTOGRAM_X_MAX PERCENTILE_BOXPLOT_Y_MAX \
        POSITIVE_RUNS_THRESHOLD POSITIVE_RUNS_CHUNK_SIZE POSITIVE_RUNS_MIN_LENGTH POSITIVE_RUNS_MAX_LENGTH POSITIVE_RUNS_PLOT_X_MAX POSITIVE_RUNS_NORMALIZATION \
        GENE_FFT_WINDOW GENE_FFT_PERIOD_MIN GENE_FFT_PERIOD_MAX GENE_FFT_RANKING_PERIODS RANDOMIZE_SEED RANDOMIZE_SEARCH_WINDOW RANDOMIZE_FALLBACK FRAG_COUNT_MIN FRAG_COUNT_MAX \
        FRAG_PLOT_MIN FRAG_PLOT_MAX HEATMAP_MIN_FRAG HEATMAP_MAX_FRAG HEATMAP_NORMALIZATION PEAK_SCORE_NORMALIZATION \
        SKIP_NRL SKIP_FRAGMENT_HEATMAP SKIP_REGION_EXTRACT SKIP_GENE_EXPRESSION \
        SKIP_TSS_EXPRESSION_QUINTILES SKIP_POSITIVE_RUNS SKIP_PEAK_SCORE_FREQUENCY COMBINE_PREREQUISITES_ONLY; do
        printf '%s\t%s\n' "$name" "${!name}"
    done
    printf 'INPUT_MODE\t%s\n' "$INPUT_MODE"
    printf 'BAM_COUNT\t%s\n' "${#BAMS[@]}"
    for index in "${!BAMS[@]}"; do
        printf 'BAM_%s\t%s\n' "$((index + 1))" "${BAMS[$index]}"
    done
    for index in "${!BAM_INPUTS[@]}"; do
        printf 'BAM_INPUT_%s\t%s\n' "$((index + 1))" "${BAM_INPUTS[$index]}"
    done
    printf 'FRAGMENT_COUNT\t%s\n' "${#FRAGMENTS[@]}"
    for index in "${!FRAGMENTS[@]}"; do
        printf 'FRAGMENT_%s\t%s\n' "$((index + 1))" "${FRAGMENTS[$index]}"
    done
    printf 'SOURCE_BAM_COUNT\t%s\n' "${#PROVENANCE_BAMS[@]}"
    for index in "${!PROVENANCE_BAMS[@]}"; do
        printf 'SOURCE_BAM_%s\t%s\n' "$((index + 1))" "${PROVENANCE_BAMS[$index]}"
    done
    printf 'SOURCE_FRAGMENT_COUNT\t%s\n' "${#PROVENANCE_FRAGMENTS[@]}"
    for index in "${!PROVENANCE_FRAGMENTS[@]}"; do
        printf 'SOURCE_FRAGMENT_%s\t%s\n' "$((index + 1))" "${PROVENANCE_FRAGMENTS[$index]}"
    done
    for index in "${!EXPRESSION_FOCUS_PROFILES[@]}"; do
        printf 'EXPRESSION_FOCUS_PROFILE_%s\t%s\n' "$((index + 1))" "${EXPRESSION_FOCUS_PROFILES[$index]}"
    done
    {
        IFS=','
        printf 'CONTIGS\t%s\n' "${CONTIGS[*]}"
    }
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

PASS_COUNT=0; FAIL_COUNT=0; SKIP_COUNT=0; FAILED_STEPS=(); FAILED_STEP_STATUS=(); FAILED_STEP_REASON=(); FAILED_STEP_LOG=()
run_step() {
    local step="$1" expected="$2"; shift 2
    local marker="$DONE_DIR/${SUPPORT_PREFIX}${PARAM_HASH}_${step}.done" log="$LOG_DIR/${SUPPORT_PREFIX}${step}.log"
    if [[ "$FORCE" -eq 0 && "$REUSE_EXISTING_OUTPUTS" -eq 1 && ( -f "$marker" || "$TRUST_EXISTING_OUTPUTS" -eq 1 ) ]]; then
        if [[ "$expected" == "-" ]] || compgen -G "$expected" >/dev/null; then
            echo "[SKIP] $step"; SKIP_COUNT=$((SKIP_COUNT + 1)); return 0
        fi
    fi
    rm -f "$marker"
    echo; echo "================================================================"
    echo "[RUN ] $step"; printf '[CMD ] '; printf '%q ' "$@"; printf '\n[LOG ] %s\n' "$log"
    echo "================================================================"
    # Run each step independently. Failures are recorded, but a failed step
    # returns success to the shell driver so all remaining steps are attempted.
    set +e
    { printf '[COMMAND] '; printf '%q ' "$@"; printf '\n\n'; "$@"; } > >(tee "$log") 2>&1
    local status=$?
    set -e
    if [[ "$status" -ne 0 ]]; then
        echo "[FAIL] $step exited with status $status" >&2
        FAIL_COUNT=$((FAIL_COUNT + 1)); FAILED_STEPS+=("$step")
        FAILED_STEP_STATUS+=("$status"); FAILED_STEP_REASON+=("command_failed"); FAILED_STEP_LOG+=("$log")
        return 0
    fi
    if [[ "$expected" != "-" ]] && ! compgen -G "$expected" >/dev/null; then
        echo "[FAIL] $step did not create expected output: $expected" >&2
        FAIL_COUNT=$((FAIL_COUNT + 1)); FAILED_STEPS+=("$step")
        FAILED_STEP_STATUS+=("0"); FAILED_STEP_REASON+=("missing_expected_output"); FAILED_STEP_LOG+=("$log")
        return 0
    fi
    touch "$marker"; echo "[PASS] $step"; PASS_COUNT=$((PASS_COUNT + 1)); return 0
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
        FAIL_COUNT=$((FAIL_COUNT + 1)); FAILED_STEPS+=("$name")
        FAILED_STEP_STATUS+=("$status"); FAILED_STEP_REASON+=("$reason"); FAILED_STEP_LOG+=("$log")
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
    # Defensive fallback for shells that reap a job without populating -p.
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
    local step="$1" expected="$2"; shift 2
    if [[ "$ASYNC_CORES" -le 1 ]]; then
        run_step "$step" "$expected" "$@"
        return 0
    fi
    local marker="$DONE_DIR/${SUPPORT_PREFIX}${PARAM_HASH}_${step}.done" log="$LOG_DIR/${SUPPORT_PREFIX}${step}.log"
    local safe_name="${step//[^A-Za-z0-9._-]/_}"
    local status_file="$ASYNC_STATUS_DIR/${SUPPORT_PREFIX}${safe_name}.$$.status"
    if [[ "$FORCE" -eq 0 && "$REUSE_EXISTING_OUTPUTS" -eq 1 && ( -f "$marker" || "$TRUST_EXISTING_OUTPUTS" -eq 1 ) ]]; then
        if [[ "$expected" == "-" ]] || compgen -G "$expected" >/dev/null; then
            echo "[SKIP] $step"; SKIP_COUNT=$((SKIP_COUNT + 1)); return 0
        fi
    fi
    rm -f "$marker"
    echo "[RUN ] $step"
    { printf '[COMMAND] '; printf '%q ' "$@"; printf '\n'; } > "$log"
    (
        worker_status=0; worker_reason="command_failed"
        if "$@" >> "$log" 2>&1; then
            worker_status=0
        else
            worker_status=$?
        fi
        if [[ "$worker_status" -eq 0 ]]; then
            if [[ "$expected" == "-" ]] || compgen -G "$expected" >/dev/null; then
                worker_reason="pass"
            else
                worker_status=2; worker_reason="missing_expected_output"
            fi
        fi
        printf '%s\t%s\n' "$worker_status" "$worker_reason" > "$status_file"
    ) &
    ASYNC_NAMES+=("$step")
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

run_step "00_cli_registry" "$SETUP_DIR/${SUPPORT_PREFIX}nucleosuite_version.txt" bash -c '
set -euo pipefail; bin="$1"; out="$2"; "$bin" --version | tee "$out"
for cmd in tracks pns coverage dyads fragment-ends dinuc-profile ww-types call-peaks fragments randomize-fragments merge-bams fragment-lengths fragment-heatmap gene-sets gene-expression tss-expression-quintiles aggregate dac nrl plot positive-runs peak-score-frequency distances region-extract resources validate-inputs mnase-suite; do "$bin" "$cmd" --help >/dev/null; done
' _ "$NUCLEOSUITE_BIN" "$SETUP_DIR/${SUPPORT_PREFIX}nucleosuite_version.txt"

run_step "00_python_dependencies" "$SETUP_DIR/${SUPPORT_PREFIX}python_dependencies.txt" "$PYTHON_BIN" - "$SETUP_DIR/${SUPPORT_PREFIX}python_dependencies.txt" <<'PY'
import importlib, sys
from pathlib import Path
mods = ["nucleosuite", "pysam", "pyBigWig", "numpy", "scipy", "matplotlib", "openpyxl"]
rows=[]
for name in mods:
    module=importlib.import_module(name); rows.append(f"{name}\t{getattr(module, '__version__', 'available')}")
Path(sys.argv[1]).write_text("module\tversion\n"+"\n".join(rows)+"\n"); print("\n".join(rows))
PY

INDEX_REPORT="$SETUP_DIR/${SUPPORT_PREFIX}input_indexes.tsv"
run_step "00_prepare_indexes" "$INDEX_REPORT" "$PYTHON_BIN" - "$INDEX_REPORT" "$FASTA" "${BAMS[@]}" <<'PY'
import sys
from pathlib import Path
import pysam
report=Path(sys.argv[1]); fasta=Path(sys.argv[2]); bams=[Path(value) for value in sys.argv[3:]]
rows=[]
for bam in bams:
    candidates=[Path(str(bam)+".bai"), bam.with_suffix(".bai"), Path(str(bam)+".csi")]
    if not any(path.exists() for path in candidates):
        pysam.index(str(bam))
    index=next((path for path in candidates if path.exists()), None)
    if index is None:
        raise SystemExit(f"Unable to create BAM index for {bam}")
    rows.append((str(bam), str(index)))
if not Path(str(fasta)+".fai").exists():
    pysam.faidx(str(fasta))
report.write_text("bam\tindex\n" + "".join(f"{bam}\t{index}\n" for bam,index in rows))
PY

# Define every observed output before setup so the manifest can be written once.
# The combined command is invoked only after 00_setup and 00_gene_sets finish.
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
PNS_SCALED_BW="$SCALED_DIR/${SAMPLE}_PNS_scaled_to_mean_nucleosome_peak_score.bw"
PNS_POS_SCALED_BW="$SCALED_DIR/${SAMPLE}_posPNS_mean_scaled.bw"
PNS_COVERAGE_SCALED_BW="$SCALED_DIR/${SAMPLE}_coverage_mean_scaled.bw"
PNS_ANALYSIS_BW="$PNS_BW"

FINE_RANGE_LABEL="${FINE_FRAG_LOWER}-${FINE_FRAG_UPPER}"
DYAD_RANGE_FOLDER="$DYAD_RANGE_DIR/$FINE_RANGE_LABEL"
DYAD_RANGE_BASE="$DYAD_RANGE_FOLDER/${SAMPLE}_${FINE_FRAG_LOWER}_${FINE_FRAG_UPPER}"
DYAD_RANGE_PREFIX="${DYAD_RANGE_BASE}_dyads_lower${FINE_FRAG_LOWER}_upper${FINE_FRAG_UPPER}"
DYAD_RANGE_BW="${DYAD_RANGE_PREFIX}_dyad.bw"
DYAD_EXACT_FOLDER="$DYAD_EXACT_DIR/$EXACT_SIZE"
DYAD_EXACT_BASE="$DYAD_EXACT_FOLDER/${SAMPLE}_${EXACT_SIZE}"
DYAD_EXACT_PREFIX="${DYAD_EXACT_BASE}_dyads_lower${EXACT_SIZE}_upper${EXACT_SIZE}"
DYAD_EXACT_BW="${DYAD_EXACT_PREFIX}_dyad.bw"
ENDS_RANGE_FOLDER="$ENDS_RANGE_DIR/$FINE_RANGE_LABEL"
ENDS_BASE="$ENDS_RANGE_FOLDER/${SAMPLE}_${FINE_FRAG_LOWER}_${FINE_FRAG_UPPER}"
ENDS_PREFIX="${ENDS_BASE}_fragment_ends_lower${FINE_FRAG_LOWER}_upper${FINE_FRAG_UPPER}"
ENDS_LEFT_BW="${ENDS_PREFIX}_fragment_left_ends.bw"
ENDS_RIGHT_BW="${ENDS_PREFIX}_fragment_right_ends.bw"
ENDS_EXACT_FOLDER="$ENDS_EXACT_DIR/$EXACT_SIZE"
ENDS_EXACT_PREFIX="$ENDS_EXACT_FOLDER/${SAMPLE}_${EXACT_SIZE}_fragment_ends_lower${EXACT_SIZE}_upper${EXACT_SIZE}"
ENDS_EXACT_LEFT_BW="${ENDS_EXACT_PREFIX}_fragment_left_ends.bw"
ENDS_EXACT_RIGHT_BW="${ENDS_EXACT_PREFIX}_fragment_right_ends.bw"

OBS_TRACK_SPEC="$COMBINED_TRACK_DIR/${SUPPORT_PREFIX}manifest.tsv"
printf 'fragment_range\toutput_prefix\ttracks\tbasic_scope\n' > "$OBS_TRACK_SPEC"
printf '%s-%s\t%s\t%s\trange\n' "$PNS_FRAG_LOWER" "$PNS_FRAG_UPPER" "$PNS_PREFIX" "$PNS_TRACK_LIST" >> "$OBS_TRACK_SPEC"
printf '%s-%s\t%s\tdyad\trange\n' "$FINE_FRAG_LOWER" "$FINE_FRAG_UPPER" "$DYAD_RANGE_PREFIX" >> "$OBS_TRACK_SPEC"
printf '%s\t%s\tdyad\trange\n' "$EXACT_SIZE" "$DYAD_EXACT_PREFIX" >> "$OBS_TRACK_SPEC"
printf '%s-%s\t%s\tfragment_ends,fragment_left_ends,fragment_right_ends\trange\n' \
    "$FINE_FRAG_LOWER" "$FINE_FRAG_UPPER" "$ENDS_PREFIX" >> "$OBS_TRACK_SPEC"
printf '%s\t%s\tfragment_left_ends,fragment_right_ends\trange\n' \
    "$EXACT_SIZE" "$ENDS_EXACT_PREFIX" >> "$OBS_TRACK_SPEC"

# MNase dinucleotide profiles retain exact 145 and 147 bp in addition to the 146-148 range.
for spec in "range:$FINE_FRAG_LOWER:$FINE_FRAG_UPPER" "exact:$DINUC_EXACT_A:$DINUC_EXACT_A" "exact:$DINUC_EXACT_B:$DINUC_EXACT_B"; do
    IFS=: read -r class lo hi <<< "$spec"
    label="$lo-$hi"
    [[ "$lo" == "$hi" ]] && label="$lo"
    dinuc_dir="$SEQUENCE_DIR/dinucleotide_profiles/${class}/${label}"
    dinuc_prefix="$dinuc_dir/${SAMPLE}_dinuc_lower${lo}_upper${hi}"
    printf '%s-%s\t%s\tdinuc_profile\trange\n' "$lo" "$hi" "$dinuc_prefix" >> "$OBS_TRACK_SPEC"
done
WW_DIR="$SEQUENCE_DIR/ww_types/ranges/$FINE_RANGE_LABEL"
TYPE_DYAD_DIR="$SEQUENCE_DIR/type_dyads/ranges/$FINE_RANGE_LABEL"
WW_PREFIX="$WW_DIR/${SAMPLE}_wwtypes_lower${FINE_FRAG_LOWER}_upper${FINE_FRAG_UPPER}"
TYPE_DYAD_PREFIX="$TYPE_DYAD_DIR/${SAMPLE}_wwtypes_lower${FINE_FRAG_LOWER}_upper${FINE_FRAG_UPPER}"
printf '%s-%s\t%s\tww_types\trange\n' "$FINE_FRAG_LOWER" "$FINE_FRAG_UPPER" "$WW_PREFIX" >> "$OBS_TRACK_SPEC"
printf '%s-%s\t%s\ttype_dyads\trange\n' "$FINE_FRAG_LOWER" "$FINE_FRAG_UPPER" "$TYPE_DYAD_PREFIX" >> "$OBS_TRACK_SPEC"
WW_SUMMARY="${WW_PREFIX}_ww_type_summary.tsv"
WW_TYPE_TRACKS=("${TYPE_DYAD_PREFIX}_type1_dyad.bw" "${TYPE_DYAD_PREFIX}_type2_dyad.bw" "${TYPE_DYAD_PREFIX}_type3_dyad.bw" "${TYPE_DYAD_PREFIX}_type4_dyad.bw")

ALL_CHROM_SIZES="$SETUP_DIR/${SUPPORT_PREFIX}analysis.chrom.sizes"
CHROM_SIZES="$SETUP_DIR/${SUPPORT_PREFIX}selected.chrom.sizes"
CONTIGS_CSV_FILE="$SETUP_DIR/${SUPPORT_PREFIX}analysis_contigs.csv"
CHROM_SIZE_INPUTS=()
if [[ -n "$ANALYSIS_CHROM_SIZES_SOURCE" ]]; then
    CHROM_SIZE_INPUTS=("$ANALYSIS_CHROM_SIZES_SOURCE")
elif [[ "$INPUT_MODE" == "bam" ]]; then
    CHROM_SIZE_INPUTS=("${BAMS[@]}")
else
    CHROM_SIZE_INPUTS=("${FRAGMENTS[@]}")
fi
run_step "00_chrom_sizes" "$CHROM_SIZES" "$PYTHON_BIN" - \
    "$ANALYSIS_CHROM_SIZES_SOURCE" "$FASTA" "$ALL_CHROM_SIZES" "$CHROM_SIZES" "$CONTIGS_CSV_FILE" \
    "${CONTIGS[*]}" "$INPUT_MODE" "${CHROM_SIZE_INPUTS[@]}" <<'PY'
import sys
from pathlib import Path
import pysam
from nucleosuite.core.bam_headers import merge_bam_reference_headers_with_aliases
from nucleosuite.core.chrom_sizes import read_chrom_sizes_source, write_chrom_sizes_table
from nucleosuite.core.fragment_inputs import IntervalFragmentSource
from nucleosuite.core.regions import expand_contig_tokens

source_path, fasta_path, all_out, selected_out, contigs_csv, tokens_text, input_mode, *inputs = sys.argv[1:]
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
Path(contigs_csv).write_text(",".join(selected) + "\n")
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
        (( FINE_FRAG_LOWER < RANDOM_LOWER )) && RANDOM_LOWER="$FINE_FRAG_LOWER"
        (( FINE_FRAG_UPPER > RANDOM_UPPER )) && RANDOM_UPPER="$FINE_FRAG_UPPER"
        (( EXACT_SIZE < RANDOM_LOWER )) && RANDOM_LOWER="$EXACT_SIZE"
        (( EXACT_SIZE > RANDOM_UPPER )) && RANDOM_UPPER="$EXACT_SIZE"
        (( DINUC_EXACT_A < RANDOM_LOWER )) && RANDOM_LOWER="$DINUC_EXACT_A"
        (( DINUC_EXACT_B > RANDOM_UPPER )) && RANDOM_UPPER="$DINUC_EXACT_B"
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

 ANALYSIS_CONTIGS_CSV="$(tr -d '\r\n' < "$CONTIGS_CSV_FILE" 2>/dev/null || true)"
CTCF_FILTERED="$SETUP_DIR/${SUPPORT_PREFIX}ctcf_compatible.bed"; CTCF_EXPANDED="$SETUP_DIR/${SUPPORT_PREFIX}ctcf_compatible_flank${CTCF_FLANK}.bed"
STATES_FILTERED="$SETUP_DIR/${SUPPORT_PREFIX}states_compatible.bed"
GENE_SET_SUMMARY=""
GENE_STATE_INTERVAL=""
if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 0 ]]; then
run_step "00_prepare_regions" "$CTCF_FILTERED" "$PYTHON_BIN" - "$CHROM_SIZES" "$CTCF_BED" "$STATES_BED" \
    "$CTCF_FILTERED" "$CTCF_EXPANDED" "$STATES_FILTERED" "$CTCF_FLANK" "$STATES_LABEL_COLUMN" <<'PY'
import sys
from pathlib import Path
(chrom_sizes, ctcf_in, states_in, ctcf_out, ctcf_expanded,
 states_out, flank_s, state_col_s)=sys.argv[1:]
flank=int(flank_s); state_col=int(state_col_s)-1
chroms={}
with open(chrom_sizes, encoding='utf-8-sig') as handle:
    for raw in handle:
        fields=raw.split()
        if len(fields) >= 2:
            chroms[fields[0]]=int(fields[1])
if not chroms:
    raise SystemExit(f'No chromosome sizes were found in {chrom_sizes}')
def resolve_chrom(name):
    if name in chroms:
        return name
    stripped=name[3:] if name.lower().startswith('chr') else name
    candidates=[stripped, 'chr'+stripped]
    if stripped.upper() in {'M','MT'}:
        candidates.extend(['M','MT','chrM','chrMT'])
    matches=[]
    for candidate in candidates:
        if candidate in chroms and candidate not in matches:
            matches.append(candidate)
    if len(matches)>1:
        raise SystemExit(
            f"Ambiguous chromosome alias {name!r} in analysis BigWig: {', '.join(matches)}"
        )
    return matches[0] if matches else None
def rows(path):
    if not path: return []
    out=[]
    with open(path, encoding='utf-8-sig') as h:
        for line_no, raw in enumerate(h,1):
            t=raw.strip()
            if not t or t.startswith(('#','track','browser')): continue
            f=t.split()
            if len(f)<3: continue
            chrom=resolve_chrom(f[0])
            if chrom is None: continue
            try: s=max(0,int(f[1])); e=min(chroms[chrom],int(f[2]))
            except ValueError: continue
            if e>s:
                f[0]=chrom
                out.append((chrom,s,e,f))
    return out
ct=[]; ex=[]
for chrom,s,e,f in rows(ctcf_in):
    strand='.'
    for index in (5,3):
        if index < len(f) and f[index] in ('+','-'):
            strand=f[index]; break
    name = f[3] if len(f)>3 and f[3] not in ('+','-','.') else f'CTCF_{len(ct)+1}'
    score='0'
    for index in (4,8):
        if index < len(f):
            try:
                float(f[index]); score=f[index]; break
            except ValueError:
                pass
    ct.append((chrom,s,e,name,score,strand)); center=(s+e)//2; a=max(0,center-flank); b=min(chroms[chrom],center+flank+1)
    if b>a: ex.append((chrom,a,b,name,score,strand))
if not ct: raise SystemExit('No compatible CTCF records remain')
for path,data in [(ctcf_out,ct),(ctcf_expanded,ex)]:
    with open(path,'w') as o:
        for row in data: o.write('\t'.join(map(str,row))+'\n')
if states_in:
    data=rows(states_in)
    if not data: raise SystemExit(f'No compatible records remain in {states_in}')
    with open(states_out,'w') as o:
        for chrom,s,e,f in data:
            # Preserve the complete ChromHMM BED record, including itemRgb in
            # column 9, while replacing clipped coordinates with the compatible
            # interval.  State labels therefore remain in their original column.
            f = list(f)
            f[1] = str(s)
            f[2] = str(e)
            o.write('\t'.join(f) + '\n')
PY
if [[ -n "$GENES_BED" && -n "$STATES_BED" ]]; then
    GENE_SET_OUTPUT_PREFIX="gene_sets"
    GENE_SET_MEMBER_ARGS=()
    if [[ "$RUN_MODE" == "randomized" ]]; then
        GENE_SET_OUTPUT_PREFIX="${SAMPLE}_gene_sets"
        GENE_SET_MEMBER_ARGS=(--prefix-member-files)
    fi
    GENE_SET_SUMMARY="$GENE_SET_DIR/${GENE_SET_OUTPUT_PREFIX}_summary.tsv"
    if [[ -s "$GENE_SET_SUMMARY" ]] && ! head -n 1 "$GENE_SET_SUMMARY" | grep -q $'\tfinal_tss_interval'; then
        echo "[INFO] Rebuilding gene sets with strand-aware TSS intervals"
        rm -f "$GENE_SET_SUMMARY"
        rm -rf "$GENE_SET_DIR/final_tss" "$TSS_AGG_DIR"
        mkdir -p "$TSS_AGG_DIR"
    fi
    IFS=',' read -r -a VENN_SET_ARRAY <<< "$VENN_SETS"
    run_step "00_gene_sets" "$GENE_SET_SUMMARY" "$NUCLEOSUITE_BIN" gene-sets \
        "${BLACKLIST_ARGS[@]}" --genes-bed "$GENES_BED" --states-bed "$STATES_BED" --config "$GENE_SET_CONFIG" \
        --state-label-column "$STATES_LABEL_COLUMN" --chrom-sizes "$CHROM_SIZES" \
        --output-dir "$GENE_SET_DIR" --output-prefix "$GENE_SET_OUTPUT_PREFIX" \
        --leftover-set-name leftover_genes --venn-sets "${VENN_SET_ARRAY[@]}" \
        "${GENE_SET_MEMBER_ARGS[@]}" --interval-format "$INTERVAL_FORMAT"
    GENE_STATE_INTERVAL="$GENE_SET_DIR/${GENE_SET_OUTPUT_PREFIX}_final_states.bed"
    if [[ "$INTERVAL_FORMAT" == "bigbed" ]]; then
        GENE_STATE_INTERVAL="$GENE_SET_DIR/${GENE_SET_OUTPUT_PREFIX}_final_states.bb"
    fi
    [[ -f "$GENE_STATE_INTERVAL" ]] || fatal "Gene-category interval was not created: $GENE_STATE_INTERVAL"
fi
fi

OBS_TRACK_REPORT="$COMBINED_TRACK_DIR/${SUPPORT_PREFIX}completion_report.tsv"
run_step "01_combined_tracks" "$OBS_TRACK_REPORT" "$NUCLEOSUITE_BIN" tracks \
    "${ANALYSIS_INPUT_ARGS[@]}" "${BLACKLIST_ARGS[@]}" --fasta "$FASTA" --chrom-sizes "$CHROM_SIZES" -c "${CONTIGS[@]}" \
    --spec-file "$OBS_TRACK_SPEC" --max-duplicates "$ACTIVE_MAX_DUPLICATES" \
    --max-per-coordinate "$MAX_PER_COORDINATE" --dedup-scope "$DEDUP_SCOPE" \
    --even-dyad "$EVEN_DYAD" --pns-mode-length "$PNS_MODE_LENGTH" --bigbed-score-scale "$BIGBED_SCORE_SCALE" \
    --pns-smooth-window "$PNS_SMOOTH_WINDOW" --pns-smooth-order "$PNS_SMOOTH_ORDER" \
    --pns-max-neg-run "$PNS_MAX_NEG_RUN" --interval-format "$INTERVAL_FORMAT" --output-format bigwig \
    --report "$OBS_TRACK_REPORT"


# Sequence profiles, WW/SS types and type-specific dyads were written by
# 01_combined_tracks. Downstream analyses use those paths directly.

# Peak calls are outputs of 01_combined_tracks and are not recalled separately.
PNS_CALL_PREFIX="$PNS_PREFIX"
PNS_CALL_NUC="${PNS_PREFIX}_nucleosome_regions.${INTERVAL_EXT}"
PNS_CALL_BRK="${PNS_PREFIX}_breakpoint_peaks.${INTERVAL_EXT}"
[[ -s "$PNS_CALL_NUC" ]] || die "PNS nucleosome peak file was not created: $PNS_CALL_NUC"
[[ -s "$PNS_CALL_BRK" ]] || die "PNS breakpoint peak file was not created: $PNS_CALL_BRK"

if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 0 ]]; then

# Scaling is deliberately post-combine so all selected chromosomes share one reference.
mkdir -p "$SCALED_DIR"
run_step "01_scale_coverage" "$PNS_COVERAGE_SCALED_BW" "$NUCLEOSUITE_BIN" mean-scale \
    --bigwig "$PNS_COVERAGE_BW" --scale 100 --output "$PNS_COVERAGE_SCALED_BW"
run_step "01_scale_pospns" "$PNS_POS_SCALED_BW" "$NUCLEOSUITE_BIN" mean-scale \
    --bigwig "$PNS_POS_BW" --scale 100 --output "$PNS_POS_SCALED_BW"
run_step "01_scale_pns_peak_mean" "$PNS_SCALED_BW" "$NUCLEOSUITE_BIN" mean-scale \
    --bigwig "$PNS_BW" --regions "$PNS_CALL_NUC" --score-column 5 --scale 100 --output "$PNS_SCALED_BW"
PNS_ANALYSIS_BW="$PNS_SCALED_BW"

DISTANCE_COMPARE_TICK_ARGS=(
    --score-z-limit "$SCORE_Z_LIMIT"
    --histogram-x-max "$DISTANCE_HISTOGRAM_X_MAX"
    --percentile-boxplot-y-max "$PERCENTILE_BOXPLOT_Y_MAX"
)
DISTANCE_NRL_TICK_ARGS=()
DISTANCE_STATE_TICK_ARGS=()
if [[ -n "$DISTANCE_X_MAJOR_TICK" ]]; then
    DISTANCE_COMPARE_TICK_ARGS+=(--distance-x-major-tick "$DISTANCE_X_MAJOR_TICK")
    DISTANCE_NRL_TICK_ARGS+=(--x-major-tick "$DISTANCE_X_MAJOR_TICK")
    DISTANCE_STATE_TICK_ARGS+=(--state-overlay-x-major-tick "$DISTANCE_X_MAJOR_TICK")
fi
if [[ -n "$DISTANCE_X_MINOR_TICK" ]]; then
    DISTANCE_COMPARE_TICK_ARGS+=(--distance-x-minor-tick "$DISTANCE_X_MINOR_TICK")
    DISTANCE_NRL_TICK_ARGS+=(--x-minor-tick "$DISTANCE_X_MINOR_TICK")
    DISTANCE_STATE_TICK_ARGS+=(--state-overlay-x-minor-tick "$DISTANCE_X_MINOR_TICK")
fi

# 02_dac: DAC is calculated only from the ranged 146-148 bp dyad track.
run_dac_scope() {
    local step="$1" track="$2" output_dir="$3" prefix="$4" regions="${5:-}" state="${6:-Combined chromosomes}"
    mkdir -p "$output_dir"
    if [[ -z "$regions" ]]; then
        queue_memory_step "$step" "$output_dir/${prefix}*_DAC_*.tsv" "$NUCLEOSUITE_BIN" dac \
            "${BLACKLIST_ARGS[@]}" --bigwig "$track" --chrom-sizes "$CHROM_SIZES" --scope combined_chromosomes \
            --window-size "$DAC_WINDOW_SIZE" --dmax "$DAC_DMAX" --algorithm "$DAC_ALGORITHM" \
            --out-prefix "$prefix" --output-dir "$output_dir" --progress-every 100
    else
        queue_memory_step "$step" "$output_dir/${prefix}*_DAC_*.tsv" "$NUCLEOSUITE_BIN" dac \
            "${BLACKLIST_ARGS[@]}" --bigwig "$track" --regions-bed "$regions" --state-column 4 --state-name "$state" \
            --dmax "$DAC_DMAX" --algorithm "$DAC_ALGORITHM" --out-prefix "$prefix" \
            --output-dir "$output_dir" --progress-every 100
    fi
}
DAC_RANGE_DIR="$DAC_DIR/dyads/ranges/$FINE_RANGE_LABEL"
run_dac_scope "02_dac_dyad_range_${FINE_FRAG_LOWER}_${FINE_FRAG_UPPER}_combined_chromosomes" "$DYAD_RANGE_BW" \
    "$DAC_RANGE_DIR/combined_chromosomes" "${SAMPLE}_dyad_${FINE_FRAG_LOWER}_${FINE_FRAG_UPPER}_combined_chromosomes"
if [[ -n "$GENE_STATE_INTERVAL" && -f "$GENE_STATE_INTERVAL" ]]; then
    run_dac_scope "02_dac_dyad_range_${FINE_FRAG_LOWER}_${FINE_FRAG_UPPER}_gene_sets" "$DYAD_RANGE_BW" \
        "$DAC_RANGE_DIR/gene_sets" "${SAMPLE}_dyad_${FINE_FRAG_LOWER}_${FINE_FRAG_UPPER}_gene_sets" \
        "$GENE_STATE_INTERVAL" gene_sets
fi
wait_queued_steps

DAC_VALIDATION_DIR="$DAC_DIR/validation"
mkdir -p "$DAC_VALIDATION_DIR"
run_step "02_verify_dac_columns" "$DAC_VALIDATION_DIR/${SUPPORT_PREFIX}DAC_COLUMN_VALIDATION.tsv" "$PYTHON_BIN" - \
    "$DAC_DIR" "$DAC_VALIDATION_DIR/${SUPPORT_PREFIX}DAC_COLUMN_VALIDATION.tsv" <<'PYDAC'
import csv, sys
from pathlib import Path
folder,out=map(Path,sys.argv[1:]); required={'DAC Value','DAC Value Percent','Raw DAC Value','Opportunities'}; rows=[]
for path in sorted(folder.rglob('*DAC*.tsv')):
    if path.name.endswith(('_summary.tsv', 'DAC_COLUMN_VALIDATION.tsv')): continue
    with path.open() as h: header=set(next(csv.reader(h,delimiter='\t')))
    missing=sorted(required-header); rows.append((path,'PASS' if not missing else 'FAIL',','.join(missing)))
if not rows: raise SystemExit('No DAC result TSVs found')
out.write_text('file\tstatus\tmissing_columns\n'+''.join(f'{p}\t{s}\t{m}\n' for p,s,m in rows))
if any(s=='FAIL' for _,s,_ in rows): raise SystemExit('One or more DAC files lack required columns')
PYDAC

# 04_nrl: three periodicity windows from the ranged-dyad DAC curves.
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
        --skip-first-peaks "$skip_first" "${DISTANCE_NRL_TICK_ARGS[@]}" \
        --output-prefix "$output_prefix" --title "$input_stem"
}
if [[ "$SKIP_NRL" -eq 0 ]]; then
    while IFS= read -r -d '' input_tsv; do
        [[ "$input_tsv" == *_summary.tsv ]] && continue
        [[ "$(basename "$input_tsv")" == *DAC_COLUMN_VALIDATION.tsv ]] && continue
        relative="${input_tsv#"$DAC_RANGE_DIR"/}"; parent="$(dirname "$relative")"
        outdir="$NRL_DIR/from_dac/dyads/ranges/$FINE_RANGE_LABEL/$parent"
        run_nrl_analysis "$input_tsv" "$outdir" "nrl_${NRL_MIN_DISTANCE}_${NRL_MAX_DISTANCE}" \
            "$NRL_MIN_DISTANCE" "$NRL_MAX_DISTANCE" "$NRL_PEAK_RESOLUTION" 1
        run_nrl_analysis "$input_tsv" "$outdir" "periodicity_${SHORT_PERIODICITY_MIN}_${SHORT_PERIODICITY_MAX}" \
            "$SHORT_PERIODICITY_MIN" "$SHORT_PERIODICITY_MAX" 1 0
        run_nrl_analysis "$input_tsv" "$outdir" "periodicity_${INTERMEDIATE_PERIODICITY_MIN}_${INTERMEDIATE_PERIODICITY_MAX}" \
            "$INTERMEDIATE_PERIODICITY_MIN" "$INTERMEDIATE_PERIODICITY_MAX" "$INTERMEDIATE_PERIODICITY_RESOLUTION" 0
    done < <(find "$DAC_RANGE_DIR" -type f -name "*DAC*.tsv" -print0 | sort -z)
fi
wait_queued_steps

# 05_ctcf_aggregation: one directory for each input track.
aggregate_track() {
    local step="$1" track="$2" output_dir="$3" prefix="$4" colour="$5" ylabel="$6"
    mkdir -p "$output_dir"
    queue_step "$step" "$output_dir/${prefix}_win*_heatmap.${PLOT_EXT}" "$NUCLEOSUITE_BIN" aggregate \
        "${BLACKLIST_ARGS[@]}" --bigwig "$track" --region-bed "$CTCF_FILTERED" --output-dir "$output_dir" --output-prefix "$prefix" \
        --window-half "$AGGREGATE_WINDOW_HALF" --strand-col 6 --missing-strand error --zero-thresh 0 \
        --max-score 1000000000000 --nan-to-zero --sort-mode mean_absolute \
        --colorbar-label "$colour" --mean-ylabel "$ylabel"
}
aggregate_track "05_ctcf_pns" "$PNS_ANALYSIS_BW" "$CTCF_AGG_DIR/pns" "${SAMPLE}_CTCF_PNS" PNS "Mean PNS"
aggregate_track "05_ctcf_dyad_range" "$DYAD_RANGE_BW" "$CTCF_AGG_DIR/dyads/ranges/$FINE_RANGE_LABEL" \
    "${SAMPLE}_CTCF_dyad_${FINE_FRAG_LOWER}_${FINE_FRAG_UPPER}" Dyad "Mean dyad signal"
aggregate_track "05_ctcf_dyad_exact" "$DYAD_EXACT_BW" "$CTCF_AGG_DIR/dyads/exact/$EXACT_SIZE" \
    "${SAMPLE}_CTCF_dyad_${EXACT_SIZE}" Dyad "Mean ${EXACT_SIZE} bp dyad signal"
for idx in 1 2 3 4; do
    aggregate_track "05_ctcf_type_dyad_${idx}" "${WW_TYPE_TRACKS[$((idx-1))]}" \
        "$CTCF_AGG_DIR/type_dyads/ranges/$FINE_RANGE_LABEL/type${idx}" \
        "${SAMPLE}_CTCF_type${idx}_dyad" Dyad "Mean type ${idx} dyad signal"
done

wait_queued_steps

# 06_tss_aggregation: organised first by input track, then gene set.
tss_aggregate_track() {
    local label="$1" track="$2" track_dir="$3" ylabel="$4"
    [[ -s "$GENE_SET_SUMMARY" ]] || return 0
    mkdir -p "$track_dir"
    local set_name tss_interval prefix profile mean_png set_dir
    local -a profile_specs=() sparse_args=(--zero-thresh 0 --nan-to-zero)
    if [[ "$label" == dyad_* || "$label" == type* ]]; then sparse_args+=(--max-score 1000000000000); fi
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
    done < <("$PYTHON_BIN" - "$GENE_SET_SUMMARY" <<'PYTSS'
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
        combined_dir="$track_dir/combined_plots"; mkdir -p "$combined_dir"
        combined_prefix="${SAMPLE}_${label}_TSS_gene_sets_combined"
        run_step "06_tss_overlay_${label}" "$combined_dir/${combined_prefix}.${PLOT_EXT}" "$PYTHON_BIN" - \
            "$combined_dir/${combined_prefix}.tsv" "$combined_dir/${combined_prefix}.${PLOT_EXT}" \
            "${SAMPLE}: ${label} at gene-set TSS" "$ylabel" "${profile_specs[@]}" <<'PYTSSPLOT'
import sys
from nucleosuite.profile_plots import plot_profile_overlay
output_tsv,output_png,title,ylabel,*specs=sys.argv[1:]
plot_profile_overlay([tuple(s.split('=',1)) for s in specs],output_tsv,output_png,xlabel='Position relative to TSS (bp)',ylabel=ylabel,title=title)
PYTSSPLOT
    fi
}
tss_aggregate_track pns "$PNS_ANALYSIS_BW" "$TSS_AGG_DIR/pns" "Mean PNS"
tss_aggregate_track "dyad_range_${FINE_RANGE_LABEL}" "$DYAD_RANGE_BW" "$TSS_AGG_DIR/dyads/ranges/$FINE_RANGE_LABEL" "Mean dyad signal"
tss_aggregate_track "dyad_exact_${EXACT_SIZE}" "$DYAD_EXACT_BW" "$TSS_AGG_DIR/dyads/exact/$EXACT_SIZE" "Mean ${EXACT_SIZE} bp dyad signal"
for idx in 1 2 3 4; do
    tss_aggregate_track "type${idx}_dyad" "${WW_TYPE_TRACKS[$((idx-1))]}" \
        "$TSS_AGG_DIR/type_dyads/ranges/$FINE_RANGE_LABEL/type${idx}" "Mean type ${idx} dyad signal"
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
source_dir="$DIST_DIR/pns_peaks"
mkdir -p "$source_dir/combined_chromosomes"
adjacent_prefix="$source_dir/combined_chromosomes/${SAMPLE}_PNS_peak_distances_adjacent"
queue_memory_step "07_distances_pns_adjacent" "${adjacent_prefix}*.tsv" "$NUCLEOSUITE_BIN" distances "$PNS_CALL_NUC" \
    "${BLACKLIST_ARGS[@]}" --position-column 7 --score-column 5 --score-percentile 0 --min-distance 1 --max-distance "$DISTANCE_ADJACENT_MAX" \
    --max-order 1 --scope combined_chromosomes --write-filtered-bed --interval-format "$INTERVAL_FORMAT" \
    --interval-chrom-sizes "$CHROM_SIZES" --output-prefix "$adjacent_prefix"
long_prefix="$source_dir/combined_chromosomes/${SAMPLE}_PNS_peak_distances_orders1-${DISTANCE_LONG_MAX_ORDER}"
queue_memory_step "07_distances_pns_nrl" "${long_prefix}*.tsv" "$NUCLEOSUITE_BIN" distances "$PNS_CALL_NUC" \
    "${BLACKLIST_ARGS[@]}" --position-column 7 --score-column 5 --score-percentile 0 --min-distance 1 --max-distance "$DISTANCE_LONG_MAX" \
    --max-order "$DISTANCE_LONG_MAX_ORDER" --scope combined_chromosomes --regression-scope combined \
    --write-filtered-bed --interval-format "$INTERVAL_FORMAT" --interval-chrom-sizes "$CHROM_SIZES" --output-prefix "$long_prefix"
if [[ -n "$STATES_BED" ]]; then
    mkdir -p "$source_dir/chromhmm_states"
    state_prefix="$source_dir/chromhmm_states/${SAMPLE}_PNS_ChromHMM_peak_distances"
    queue_memory_step "07_state_distances_pns" "${state_prefix}_scorepct0_state_relative_percent.${PLOT_EXT}" \
        "$NUCLEOSUITE_BIN" distances "$PNS_CALL_NUC" --position-column 7 --score-column 5 --score-percentile 0 \
        "${BLACKLIST_ARGS[@]}" --min-distance 1 --max-distance "$STATE_DISTANCE_MAX" --max-order 1 --scope combined_chromosomes \
        --state-bed "$STATES_FILTERED" --state-label-column "$STATES_LABEL_COLUMN" --state-color-column 9 \
        --state-overlay-plot --state-overlay-smooth-window "$STATE_DISTANCE_SMOOTH_WINDOW" \
        --state-overlay-smooth-polyorder "$STATE_DISTANCE_SMOOTH_ORDER" "${DISTANCE_STATE_TICK_ARGS[@]}" \
        --state-overlay-title "${SAMPLE} PNS: adjacent peak distances by ChromHMM state" --output-prefix "$state_prefix"
fi
wait_queued_steps

# 08_region_extract: CTCF extractions are grouped by signal track.
if [[ "$SKIP_REGION_EXTRACT" -eq 0 ]]; then
    PNS_REGION_DIR="$REGION_DIR/ctcf/pns"; mkdir -p "$PNS_REGION_DIR"
    PNS_REGION_PREFIX="$PNS_REGION_DIR/${SAMPLE}_CTCF_PNS"
    queue_step "08_region_extract_pns" "${PNS_REGION_PREFIX}_pns_signal.tsv" "$NUCLEOSUITE_BIN" region-extract \
        "${BLACKLIST_ARGS[@]}" --bed "$CTCF_EXPANDED" --coverage-bw "$PNS_COVERAGE_SCALED_BW" --pns-bw "$PNS_ANALYSIS_BW" \
        --nucleosome-peaks "$PNS_CALL_NUC" --breakpoint-peaks "$PNS_CALL_BRK" --peak-flank-bp "$REGION_PEAK_FLANK" \
        --peak-center-column 7 --peak-score-column 5 --out-prefix "$PNS_REGION_PREFIX" --chrom-mode auto \
        --missing-chrom error --progress-every 100 --overwrite
fi
wait_queued_steps

# 09_fragment_lengths: organised by the regions used for counting.
WHOLE_FRAG_DIR="$FRAG_DIR/combined_chromosomes"; mkdir -p "$WHOLE_FRAG_DIR"
WHOLE_FRAG_TSV="$WHOLE_FRAG_DIR/${SAMPLE}_combined_chromosomes_fragment_lengths.tsv"
queue_step "09_fragment_lengths_combined_chromosomes" "$WHOLE_FRAG_TSV" "$NUCLEOSUITE_BIN" fragment-lengths \
    "${ANALYSIS_INPUT_ARGS[@]}" "${BLACKLIST_ARGS[@]}" --contigs "${ANALYSIS_CONTIGS_CSV:-all}" --min-length "$FRAG_COUNT_MIN" \
    --max-length "$FRAG_COUNT_MAX" --output "$WHOLE_FRAG_TSV" \
    --plot "$WHOLE_FRAG_DIR/${SAMPLE}_combined_chromosomes_fragment_lengths.${PLOT_EXT}" --plot-min "$FRAG_PLOT_MIN" --plot-max "$FRAG_PLOT_MAX"
STATE_FRAG_TSV=""
if [[ -n "$STATES_BED" ]]; then
    STATE_FRAG_DIR="$FRAG_DIR/chromhmm_states"; mkdir -p "$STATE_FRAG_DIR"
    STATE_FRAG_TSV="$STATE_FRAG_DIR/${SAMPLE}_states_fragment_lengths.tsv"
    queue_step "09_fragment_lengths_chromhmm_states" "$STATE_FRAG_TSV" "$NUCLEOSUITE_BIN" fragment-lengths \
        "${ANALYSIS_INPUT_ARGS[@]}" "${BLACKLIST_ARGS[@]}" --bed "$STATES_FILTERED" --bed-label-column "$STATES_LABEL_COLUMN" \
        --overlap-policy all --contigs "${ANALYSIS_CONTIGS_CSV:-all}" --min-length "$FRAG_COUNT_MIN" \
        --max-length "$FRAG_COUNT_MAX" --output "$STATE_FRAG_TSV" --separate-files \
        --output-dir "$STATE_FRAG_DIR/state_profiles" --plot "$STATE_FRAG_DIR/${SAMPLE}_states_fragment_lengths.${PLOT_EXT}" \
        --plot-min "$FRAG_PLOT_MIN" --plot-max "$FRAG_PLOT_MAX"
fi
wait_queued_steps

# 10_fragment_heatmaps: combined combined-chromosome and region-specific profiles.
if [[ "$SKIP_FRAGMENT_HEATMAP" -eq 0 ]]; then
    MISSING_HEATMAP_DEPENDENCIES=()
    [[ -s "$WHOLE_FRAG_TSV" ]] || MISSING_HEATMAP_DEPENDENCIES+=("$WHOLE_FRAG_TSV")
    if [[ -n "$STATE_FRAG_TSV" && ! -s "$STATE_FRAG_TSV" ]]; then
        MISSING_HEATMAP_DEPENDENCIES+=("$STATE_FRAG_TSV")
    fi
    if [[ "${#MISSING_HEATMAP_DEPENDENCIES[@]}" -gt 0 ]]; then
        echo "[SKIP] 10_fragment_heatmap_combined: missing fragment-length dependency: ${MISSING_HEATMAP_DEPENDENCIES[*]}"
        SKIP_COUNT=$((SKIP_COUNT + 1))
    else
        HEATMAP_COMBINED_DIR="$HEATMAP_DIR/combined"; mkdir -p "$HEATMAP_COMBINED_DIR"
        HEATMAP_PREFIX="$HEATMAP_COMBINED_DIR/${SAMPLE}_fragment_lengths"
        HEATMAP_INPUTS=(--input "Combined_chromosomes=$WHOLE_FRAG_TSV")
        [[ -z "$STATE_FRAG_TSV" ]] || HEATMAP_INPUTS+=(--input "$STATE_FRAG_TSV")
        run_step "10_fragment_heatmap_combined" "${HEATMAP_PREFIX}_fragmin*_heatmap.${PLOT_EXT}" "$NUCLEOSUITE_BIN" fragment-heatmap \
            "${HEATMAP_INPUTS[@]}" --out-prefix "$HEATMAP_PREFIX" --min-frag "$HEATMAP_MIN_FRAG" \
            --max-frag "$HEATMAP_MAX_FRAG" --normalization "$HEATMAP_NORMALIZATION" \
            --title "${SAMPLE}: MNase fragment-length profiles"
    fi
fi

# 11_gene_expression: PNS-only gene-expression analysis.
if [[ -n "$EXPRESSION" && "$SKIP_GENE_EXPRESSION" -eq 0 ]]; then
    FOCUS_PROFILE_ARGS=(); for profile in "${EXPRESSION_FOCUS_PROFILES[@]}"; do FOCUS_PROFILE_ARGS+=(--focus-profile "$profile"); done
    run_gene_expression() {
        local label="$1" signal_type="$2" peaks="$3" signal="$4" output_dir="$5"
        mkdir -p "$output_dir"
        local prefix="$output_dir/${SAMPLE}_${label}_gene_expression"
        queue_memory_step "11_gene_expression_${label,,}" "${prefix}_analysis*_metadata.tsv" "$NUCLEOSUITE_BIN" gene-expression \
            "${BLACKLIST_ARGS[@]}" --expression "$EXPRESSION" --genes-bed "$GENES_BED" --peaks "${SAMPLE}=${peaks}" \
            --signal "${SAMPLE}=${signal}" --signal-type "$signal_type" --analysis all --output-prefix "$prefix" \
            --expression-gene-column "$EXPRESSION_GENE_COLUMN" --expression-name-column "$EXPRESSION_NAME_COLUMN" \
            --expression-profile-column "$EXPRESSION_PROFILE_COLUMN" --expression-value-column "$EXPRESSION_VALUE_COLUMN" \
            --fft-window "$GENE_FFT_WINDOW" --fft-period-min "$GENE_FFT_PERIOD_MIN" --fft-period-max "$GENE_FFT_PERIOD_MAX" \
            --fft-ranking-periods "$GENE_FFT_RANKING_PERIODS" "${FOCUS_PROFILE_ARGS[@]}"
    }
    run_gene_expression PNS pns "$PNS_CALL_NUC" "$PNS_ANALYSIS_BW" "$GENE_EXPRESSION_DIR/pns"
fi
wait_queued_steps

# 12_positive_runs: the active suite input (observed or randomized-only).
if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 0 && "$SKIP_POSITIVE_RUNS" -eq 0 ]]; then
    PNS_POSITIVE_DIR="$POSITIVE_RUNS_DIR/pns"
    mkdir -p "$PNS_POSITIVE_DIR"
    PNS_POSITIVE_PREFIX="$PNS_POSITIVE_DIR/${SAMPLE}_PNS_positive_runs"
    queue_step "12_positive_runs_pns" "${PNS_POSITIVE_PREFIX}_threshold*_summary.tsv" "$NUCLEOSUITE_BIN" positive-runs \
        "${BLACKLIST_ARGS[@]}" --bigwig "$PNS_ANALYSIS_BW" --output-prefix "$PNS_POSITIVE_PREFIX" --contigs "${CONTIGS[@]}" \
        --threshold "$POSITIVE_RUNS_THRESHOLD" --chunk-size "$POSITIVE_RUNS_CHUNK_SIZE" \
        --min-run-length "$POSITIVE_RUNS_MIN_LENGTH" --max-run-length "$POSITIVE_RUNS_MAX_LENGTH" \
        --plot-x-max "$POSITIVE_RUNS_PLOT_X_MAX" --normalization "$POSITIVE_RUNS_NORMALIZATION" \
        --title "${SAMPLE}: PNS positive run lengths"
fi
wait_queued_steps

# 13_peak_analysis: PNS peak-score distributions.
run_peak_score_frequency() {
    local step="$1" output_dir="$2" label="$3" peaks="$4" title="$5"
    mkdir -p "$output_dir"
    local prefix output
    prefix="$output_dir/${SAMPLE}_${label}"
    output="${prefix}_bins*_score_frequency.tsv"
    local -a peak_args=(--peaks "${RUN_MODE}=$peaks")
    queue_memory_step "$step" "$output" "$NUCLEOSUITE_BIN" peak-score-frequency "${peak_args[@]}" \
        "${BLACKLIST_ARGS[@]}" --output-prefix "$prefix" --score-column 5 --integer-bins \
        --normalization "$PEAK_SCORE_NORMALIZATION" --title "$title"
}
if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 0 && "$SKIP_PEAK_SCORE_FREQUENCY" -eq 0 ]]; then
    run_peak_score_frequency "13_peak_scores_pns_nucleosome" "$PEAK_SCORE_DIR/pns" PNS_nucleosome "$PNS_CALL_NUC" "${SAMPLE}: PNS nucleosome-region scores"
    run_peak_score_frequency "13_peak_scores_pns_breakpoint" "$PEAK_SCORE_DIR/pns" PNS_breakpoint "$PNS_CALL_BRK" "${SAMPLE}: PNS breakpoint-peak scores"
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

REPORT="$OUTDIR/${SUPPORT_PREFIX}NUCLEOSUITE_MNASE_SUITE_REPORT.tsv"
{
    INPUT_COUNT=${#BAMS[@]}
    [[ "$INPUT_MODE" == "bam" ]] || INPUT_COUNT=${#FRAGMENTS[@]}
    echo -e "metric\tvalue"; echo -e "sample\t$SAMPLE"; echo -e "input_mode\t$INPUT_MODE"; echo -e "input_count\t$INPUT_COUNT"; echo -e "run_mode\t$RUN_MODE"; echo -e "parameter_hash\t$PARAM_HASH"
    echo -e "blacklist_bed\t$BLACKLIST_BED"
    echo -e "expression_table\t${EXPRESSION:-}"; echo -e "gene_expression_signals\tpns"
    if [[ "$COMBINE_PREREQUISITES_ONLY" -eq 1 ]]; then echo -e "execution_scope\tcombine_prerequisites_only"; else echo -e "execution_scope\tcombined_chromosomes_analysis"; fi
    echo -e "passed_steps\t$PASS_COUNT"
    echo -e "failed_steps\t$FAIL_COUNT"; echo -e "skipped_completed_steps\t$SKIP_COUNT"
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
    "schema_version": 1, "nucleosuite_version": __version__, "suite": "mnase",
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
echo; echo "NucleoSuite MNase suite complete"; echo "Passed: $PASS_COUNT  Failed: $FAIL_COUNT  Skipped: $SKIP_COUNT"; echo "Report: $REPORT"; echo "Logs: $LOG_DIR"
if [[ "$FAIL_COUNT" -gt 0 ]]; then printf '  - %s\n' "${FAILED_STEPS[@]}" >&2; exit 1; fi
exit 0
