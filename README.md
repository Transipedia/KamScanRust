<!-- README.md for KamScanRust -->

# KamScanRust

**KamScanRust** Performs fast parallel statistical tests over long k-mer, contig or other features count matrices, streaming differential analysis between two conditions (Welch T-test, Wilcoxon rank-sum, or ANCOVA) — multi-threaded.
Adapted from https://github.com/Transipedia/KamScan.

---

## 🛠️ Installation

### Prerequisites

- Rust 1.70 or later ([Install Rust](https://www.rust-lang.org/tools/install))
- Git (optional)

### Steps

```bash
# Clone the repository
git clone https://github.com/Transipedia/KamScanRust.git
cd KamScanRust

# Build and install
cargo build --release
```
---

## Usage

### Basic Command

```bash
kamscan [OPTIONS] --design <DESIGN>
```

## Options

-  -t, --tsv <TSV>
          Input counts file, space- or tab-separated. Use '-' or omit to read from stdin

-  -d, --design <DESIGN>
          Design file: sample_name  condition  [covariate]

-  -o, --output <OUTPUT>
          Output file  [default: stdout]

-  -p, --pseudo <PSEUDO>
          Pseudocount added before log2 transformation  [default: 1.0]

-  -m, --min-count <MIN_COUNT>
          Skip features where ALL relevant samples have count <= this value  [default: 0.0]

-  -s, --max-pvalue <MAX_PVALUE>
          Only report features with p-value (or adj p-value when --bh) <= this threshold [default: 0.05]

-  --test <TEST>
          Statistical test to use  [default: ttest]
          Possible values:
          - ttest:    Welch's two-sample t-test on log2-transformed counts (default)
          - wilcoxon: Wilcoxon rank-sum test (Mann-Whitney U) on log2-transformed counts
          - ancova:   ANCOVA: remove linear covariate effect before testing group difference. Requires a numeric third column in the design file

-  --no-header <HEADER_FILE>
          Counts file has no header row; instead, read column names (one per line) from HEADER_FILE. The first name is the feature-ID column; remaining names are sample names matched against the design file in order

- -j, --threads <THREADS>
          Number of worker threads  [default: logical CPU count]

-  --chunk-factor <CHUNK_FACTOR>
          Chunk size multiplier: each chunk holds chunk_factor × 64 KiB of raw input text.  Increase for files with very long lines. [default: 16  →  1 MiB chunks]
          
-  --norm-counts <NORM_COUNTS>
          Optional file of pre-computed column sums for CPM normalization. Format: sample_name TAB count_sum  (one line per sample, any order)

-  --norm-scale <NORM_SCALE>
          Normalization factor for the CPM normalization [default: 1000000]

-  --bh
          Apply Benjamini-Hochberg FDR correction (two-pass run). Incompatible with reading the counts file from stdin

-  --output-mode <OUTPUT_MODE>
          What to output for significant features. 'stats': means, logFC, test statistic, p-value (default). 'raw': the original input line verbatim (useful for large count matrices)

          Possible values:
          - stats: Output statistics columns: means, logFC, test stat, p-value (default)
          - raw:   Output the original input line verbatim for each significant feature

-  -h, --help
          Print help (see a summary with '-h')
