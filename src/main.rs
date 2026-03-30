/// kamscan: Streaming k-mer / feature differential analysis (multi-threaded)
///
/// Usage:
///   kamscan --tsv counts.tsv --design design.tsv [OPTIONS]
///
/// Counts file format (space- or tab-separated):
///   feature_id  sample1  sample2  sample3  ...
///
/// Design file format (space- or tab-separated):
///   sample_name  condition  [covariate_value]
///   (third column required when --test ancova)
///
/// Threading model: reader → bounded channel → worker pool → ordered writer.
/// Two-pass BH correction available via --bh.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::Arc;
use std::thread;

use clap::{Parser, ValueEnum};

// ─── Constants ────────────────────────────────────────────────────────────────

const IO_BUF_BYTES:         usize = 65_536;
const CHANNEL_DEPTH_FACTOR: usize = 2;

// ─── CLI ──────────────────────────────────────────────────────────────────────

#[derive(Clone, Debug, ValueEnum)]
enum TestType {
    /// Welch's two-sample t-test on log2-transformed counts (default)
    Ttest,
    /// Wilcoxon rank-sum test (Mann-Whitney U) on log2-transformed counts
    Wilcoxon,
    /// ANCOVA: remove linear covariate effect before testing group difference.
    /// Requires a numeric third column in the design file.
    Ancova,
}

#[derive(Clone, Debug, ValueEnum)]
enum OutputMode {
    /// Output statistics columns: means, logFC, test stat, p-value (default)
    Stats,
    /// Output the original input line verbatim for each significant feature
    Raw,
}

#[derive(Parser, Debug)]
#[command(
    name = "kamscan",
    about = "Streaming differential analysis between two conditions \
             (Welch T-test, Wilcoxon rank-sum, or ANCOVA) — multi-threaded",
)]
struct Args {
    /// Input counts file, space- or tab-separated. Use '-' or omit to read from stdin.
    #[arg(short = 't', long)]
    tsv: Option<PathBuf>,

    /// Design file: sample_name  condition  [covariate]
    #[arg(short = 'd', long)]
    design: PathBuf,

    /// Output file  [default: stdout]
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,

    /// Pseudocount added before log2 transformation  [default: 1.0]
    #[arg(short = 'p', long, default_value_t = 1.0)]
    pseudo: f64,

    /// Skip features where ALL relevant samples have count <= this value  [default: 0.0]
    #[arg(short = 'm', long, default_value_t = 0.0)]
    min_count: f64,

    /// Only report features with p-value (or adj p-value when --bh) <= this threshold
    /// [default: 0.05]
    #[arg(short = 's', long, default_value_t = 0.05)]
    max_pvalue: f64,

    /// Statistical test to use  [default: ttest]
    #[arg(long, value_enum, default_value_t = TestType::Ttest)]
    test: TestType,

    /// Counts file has no header row; instead, read column names (one per line)
    /// from HEADER_FILE. The first name is the feature-ID column; remaining names
    /// are sample names matched against the design file in order.
    #[arg(long, value_name = "HEADER_FILE")]
    no_header: Option<PathBuf>,

    /// Number of worker threads  [default: logical CPU count]
    #[arg(short = 'j', long)]
    threads: Option<usize>,

    /// Chunk size multiplier: each chunk holds chunk_factor × 64 KiB of raw
    /// input text.  Increase for files with very long lines.
    /// [default: 16  →  1 MiB chunks]
    #[arg(long, default_value_t = 16)]
    chunk_factor: usize,

    /// Optional file of pre-computed column sums for CPM normalization.
    /// Format: sample_name TAB count_sum  (one line per sample, any order).
    #[arg(long)]
    norm_counts: Option<PathBuf>,

    /// Normalization factor for the CPM normalization.
    #[arg(long, default_value_t = 1_000_000.0)]
    norm_scale: f64,

    /// Apply Benjamini-Hochberg FDR correction (two-pass run).
    /// Incompatible with reading the counts file from stdin.
    #[arg(long, default_value_t = false)]
    bh: bool,

    /// What to output for significant features.
    /// 'stats': means, logFC, test statistic, p-value (default).
    /// 'raw': the original input line verbatim (useful for large count matrices).
    #[arg(long, value_enum, default_value_t = OutputMode::Stats)]
    output_mode: OutputMode,
}

// ─── Statistics: shared helpers ───────────────────────────────────────────────

#[inline]
fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

#[inline]
fn sample_variance(v: &[f64], m: f64) -> f64 {
    if v.len() < 2 { return 0.0; }
    v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() as f64 - 1.0)
}

// ─── Statistics: Welch's T-test ───────────────────────────────────────────────

fn welch_ttest(a: &[f64], b: &[f64]) -> Option<(f64, f64, f64)> {
    let na = a.len() as f64;
    let nb = b.len() as f64;
    if na < 2.0 || nb < 2.0 { return None; }

    let mean_a = mean(a);
    let mean_b = mean(b);
    let var_a  = sample_variance(a, mean_a);
    let var_b  = sample_variance(b, mean_b);
    let se2_a  = var_a / na;
    let se2_b  = var_b / nb;
    let se2    = se2_a + se2_b;

    if se2 == 0.0 { return Some((0.0, 1.0, na + nb - 2.0)); }

    let t  = (mean_a - mean_b) / se2.sqrt();
    let df = se2 * se2 / (se2_a * se2_a / (na - 1.0) + se2_b * se2_b / (nb - 1.0));
    Some((t, two_tailed_t_pvalue(t, df), df))
}

fn two_tailed_t_pvalue(t: f64, df: f64) -> f64 {
    let x = df / (df + t * t);
    (2.0 * 0.5 * regularized_incomplete_beta(df / 2.0, 0.5, x)).min(1.0)
}

fn regularized_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    if x >= 1.0 { return 1.0; }
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - regularized_incomplete_beta(b, a, 1.0 - x);
    }
    let lbeta = ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b);
    let front = (lbeta + a * x.ln() + b * (1.0 - x).ln()).exp() / a;
    front * beta_cf(a, b, x)
}

fn beta_cf(a: f64, b: f64, x: f64) -> f64 {
    const MAX_ITER: usize = 200;
    const EPS: f64 = 3e-7;
    const FPMIN: f64 = 1e-30;
    let qab = a + b; let qap = a + 1.0; let qam = a - 1.0;
    let mut c = 1.0_f64;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < FPMIN { d = FPMIN; }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=MAX_ITER {
        let m = m as f64; let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; if d.abs() < FPMIN { d = FPMIN; }
        c = 1.0 + aa / c; if c.abs() < FPMIN { c = FPMIN; }
        d = 1.0 / d; h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; if d.abs() < FPMIN { d = FPMIN; }
        c = 1.0 + aa / c; if c.abs() < FPMIN { c = FPMIN; }
        d = 1.0 / d; let delta = d * c; h *= delta;
        if (delta - 1.0).abs() < EPS { break; }
    }
    h
}

fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.99999999999980993, 676.5203681218851, -1259.1392167224028,
        771.32342877765313, -176.61502916214059, 12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
    ];
    if x < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().ln()
            - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + G + 0.5;
    for (i, &ci) in C[1..].iter().enumerate() { a += ci / (x + i as f64 + 1.0); }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

// ─── Statistics: Wilcoxon rank-sum (Mann-Whitney U) ──────────────────────────

fn wilcoxon_ranksum(a: &[f64], b: &[f64]) -> Option<(f64, f64)> {
    let n1 = a.len(); let n2 = b.len();
    if n1 < 2 || n2 < 2 { return None; }

    let mut pool: Vec<(f64, usize)> = a.iter().map(|&v| (v, 0))
        .chain(b.iter().map(|&v| (v, 1))).collect();
    pool.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));

    let n = pool.len();
    let mut ranks = vec![0.0_f64; n];
    let mut tie_correction = 0.0_f64;
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && (pool[j].0 - pool[i].0).abs()
            < f64::EPSILON * pool[i].0.abs().max(1.0) { j += 1; }
        let avg_rank = (i + j + 1) as f64 / 2.0;
        let t = (j - i) as f64;
        tie_correction += t * t * t - t;
        for k in i..j { ranks[k] = avg_rank; }
        i = j;
    }

    let w: f64 = pool.iter().zip(ranks.iter())
        .filter(|((_, g), _)| *g == 0).map(|(_, &r)| r).sum();

    let n1f = n1 as f64; let n2f = n2 as f64; let nf = (n1 + n2) as f64;
    let u      = w - n1f * (n1f + 1.0) / 2.0;
    let mean_u = n1f * n2f / 2.0;
    let var_u  = n1f * n2f / 12.0
        * ((nf + 1.0) - tie_correction / (nf * (nf - 1.0)));

    if var_u <= 0.0 { return Some((u, 1.0)); }
    let z_num = (u - mean_u).abs() - 0.5;
    let z = if z_num <= 0.0 { 0.0 } else { z_num / var_u.sqrt() };
    Some((u, (2.0 * standard_normal_upper_tail(z)).min(1.0)))
}

fn standard_normal_upper_tail(z: f64) -> f64 {
    if z < 0.0 { return 1.0 - standard_normal_upper_tail(-z); }
    0.5 * erfc(z / std::f64::consts::SQRT_2)
}

fn erfc(x: f64) -> f64 {
    if x < 0.0 { return 2.0 - erfc(-x); }
    if x == 0.0 { return 1.0; }
    if x < 3.0 { 1.0 - erf_series(x) } else { erfc_cf(x) }
}

fn erf_series(x: f64) -> f64 {
    let x2 = x * x; let mut term = x; let mut sum = x;
    for n in 1..=50usize {
        term *= -x2 / n as f64;
        sum  += term / (2 * n + 1) as f64;
        if term.abs() < 1e-15 * sum.abs() { break; }
    }
    sum * 2.0 / std::f64::consts::PI.sqrt()
}

fn erfc_cf(x: f64) -> f64 {
    let x2 = x * x; const FPMIN: f64 = 1e-300;
    let mut b = x2 + 0.5; let mut c = 1.0 / FPMIN;
    let mut d = 1.0 / b; let mut h = d;
    for i in 1..=60usize {
        let a = -(i as f64) * (i as f64 - 0.5); b += 2.0;
        d = a * d + b; if d.abs() < FPMIN { d = FPMIN; }
        c = b + a / c; if c.abs() < FPMIN { c = FPMIN; }
        d = 1.0 / d; let delta = d * c; h *= delta;
        if (delta - 1.0).abs() < 1e-15 { break; }
    }
    h * x * (-x2).exp() / std::f64::consts::PI.sqrt()
}

// ─── Statistics: ANCOVA ───────────────────────────────────────────────────────
//
// Model: y_i = β0 + β1·group_i + β2·covariate_i + ε_i
//
// group_i = 0 for cond1, 1 for cond2.
// y_i     = log2(count_i + pseudo).
//
// We solve the 3×3 normal equations (XᵀX)β = Xᵀy using Gaussian elimination
// with partial pivoting, then compute the standard error of β1 and the
// two-tailed p-value from the t-distribution with df = N−3.
//
// Returns (t_stat, p_value, df) for the group coefficient β1.
// Returns None if N < 4 (need at least 4 samples for 3 parameters + 1 df).

fn ancova(
    y:   &[f64],   // log2-transformed counts, all N samples concatenated
    grp: &[f64],   // group indicator: 0.0 or 1.0, same order as y
    cov: &[f64],   // covariate values, same order as y
) -> Option<(f64, f64, f64)> {
    let n = y.len();
    if n < 4 || grp.len() != n || cov.len() != n { return None; }

    // Build XᵀX (3×3, symmetric) and Xᵀy (3×1)
    // Columns of X: [1, group, covariate]
    let mut xtx = [[0.0_f64; 3]; 3];
    let mut xty = [0.0_f64; 3];

    for i in 0..n {
        let xi = [1.0_f64, grp[i], cov[i]];
        for r in 0..3 {
            xty[r] += xi[r] * y[i];
            for c in 0..3 {
                xtx[r][c] += xi[r] * xi[c];
            }
        }
    }

    // Augmented matrix [XᵀX | Xᵀy] for Gaussian elimination
    let mut aug = [[0.0_f64; 4]; 3];
    for r in 0..3 {
        for c in 0..3 { aug[r][c] = xtx[r][c]; }
        aug[r][3] = xty[r];
    }

    // Gaussian elimination with partial pivoting
    for col in 0..3 {
        // Find pivot
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..3 {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        if max_val < 1e-12 { return None; } // singular or near-singular
        aug.swap(col, max_row);

        let pivot = aug[col][col];
        for c in col..4 { aug[col][c] /= pivot; }

        for row in 0..3 {
            if row == col { continue; }
            let factor = aug[row][col];
            for c in col..4 { aug[row][c] -= factor * aug[col][c]; }
        }
    }

    // β = aug[*][3]
    let beta0 = aug[0][3];
    let beta1 = aug[1][3]; // group effect
    let beta2 = aug[2][3];

    // Residual sum of squares
    let rss: f64 = y.iter().zip(grp.iter()).zip(cov.iter())
        .map(|((&yi, &gi), &ci)| {
            let resid = yi - beta0 - beta1 * gi - beta2 * ci;
            resid * resid
        })
        .sum();

    let df   = (n as f64) - 3.0;
    let mse  = rss / df;
    if mse <= 0.0 { return None; }

    // Variance of β1 = MSE * (XᵀX)⁻¹[1,1]
    // We already have (XᵀX)⁻¹ in aug (reduced row echelon form gives us
    // the inverse implicitly). Rather than reconstruct it, we directly use
    // the formula: var(β1) = MSE / (SS_group_adjusted)
    // where SS_group_adjusted = Σ(gi - ĝi)² with ĝi = fitted values of
    // group from regressing group on the covariate.
    //
    // This is equivalent to the Frisch-Waugh-Lovell theorem:
    // partial out the covariate from the group indicator, then
    // var(β1) = MSE / Σ(ẽi²) where ẽi are residuals of group ~ covariate.

    // Step 1: regress group on [1, covariate] to get ẽ
    let n_f    = n as f64;
    let mean_c = cov.iter().sum::<f64>() / n_f;
    let mean_g = grp.iter().sum::<f64>() / n_f;
    let scc: f64 = cov.iter().map(|&c| (c - mean_c).powi(2)).sum();
    let scg: f64 = cov.iter().zip(grp.iter())
        .map(|(&c, &g)| (c - mean_c) * (g - mean_g)).sum();

    let (slope_gc, intercept_gc) = if scc.abs() < 1e-12 {
        (0.0, mean_g) // covariate is constant → no partial-out needed
    } else {
        let s = scg / scc;
        (s, mean_g - s * mean_c)
    };

    let ss_resid_g: f64 = grp.iter().zip(cov.iter())
        .map(|(&g, &c)| {
            let g_hat = intercept_gc + slope_gc * c;
            (g - g_hat).powi(2)
        })
        .sum();

    if ss_resid_g < 1e-12 { return None; } // group perfectly predicted by covariate

    let var_beta1 = mse / ss_resid_g;
    if var_beta1 <= 0.0 { return None; }

    let t = beta1 / var_beta1.sqrt();
    let p = two_tailed_t_pvalue(t, df);
    Some((t, p, df))
}

// ─── Benjamini-Hochberg correction ────────────────────────────────────────────

fn benjamini_hochberg(
    pvals: &[(String, f64)],
    q_threshold: f64,
) -> HashMap<String, f64> {
    let m = pvals.len();
    if m == 0 { return HashMap::new(); }

    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| pvals[a].1.partial_cmp(&pvals[b].1)
        .unwrap_or(std::cmp::Ordering::Equal));

    let mut adj = vec![0.0_f64; m];
    let mf = m as f64;
    for (rank, &idx) in order.iter().enumerate() {
        adj[rank] = pvals[idx].1 * mf / (rank + 1) as f64;
    }
    for i in (0..m - 1).rev() {
        if adj[i] > adj[i + 1] { adj[i] = adj[i + 1]; }
    }

    let mut result = HashMap::new();
    for (rank, &idx) in order.iter().enumerate() {
        let q = adj[rank].min(1.0);
        if q <= q_threshold { result.insert(pvals[idx].0.clone(), q); }
    }
    result
}

// ─── Field splitting ──────────────────────────────────────────────────────────

fn split_fields(s: &str) -> Vec<&str> {
    let tab_fields: Vec<&str> = s.split('\t').collect();
    let spc_fields: Vec<&str> = s.split_whitespace().collect();
    if spc_fields.len() > tab_fields.len() { spc_fields } else { tab_fields }
}

fn split_design_line(s: &str) -> Vec<&str> {
    if s.contains('\t') {
        s.split('\t').map(|f| f.trim()).collect()
    } else {
        s.split_whitespace().collect()
    }
}

// ─── Design file parsing ──────────────────────────────────────────────────────

struct Design {
    sample_to_cond:      HashMap<String, String>,
    /// Covariate value per sample (None if not present in design file)
    sample_to_covariate: HashMap<String, f64>,
    samples_ordered:     Vec<String>,
    cond1: String,
    cond2: String,
}

fn parse_design(path: &PathBuf) -> io::Result<Design> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut sample_to_cond:      HashMap<String, String> = HashMap::new();
    let mut sample_to_covariate: HashMap<String, f64>    = HashMap::new();
    let mut samples_ordered:     Vec<String>             = Vec::new();
    let mut conditions:          Vec<String>             = Vec::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?.trim().to_string();
        if line.is_empty() { continue; }
        let parts = split_design_line(&line);
        if parts.len() < 2 {
            eprintln!("Warning: design line {} has < 2 fields, skipping.", lineno + 1);
            continue;
        }
        let sample = parts[0].trim().to_string();
        let cond   = parts[1].trim().to_string();

        // Optional third column: covariate
        if parts.len() >= 3 {
            let cov_str = parts[2].trim();
            match cov_str.parse::<f64>() {
                Ok(v)  => { sample_to_covariate.insert(sample.clone(), v); }
                Err(_) => {
                    return Err(io::Error::new(io::ErrorKind::InvalidData,
                        format!("design line {}: cannot parse covariate '{}' as a number",
                            lineno + 1, cov_str)));
                }
            }
        }

        if sample_to_cond.contains_key(&sample) {
            eprintln!("Warning: duplicate sample '{}' at line {}; overwriting.", sample, lineno + 1);
        } else {
            samples_ordered.push(sample.clone());
        }
        sample_to_cond.insert(sample, cond.clone());
        if !conditions.contains(&cond) { conditions.push(cond.clone()); }
        if conditions.len() > 2 {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Design file has more than 2 conditions: {:?}", conditions)));
        }
    }
    if conditions.len() < 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("Design file must contain exactly 2 conditions; found: {:?}", conditions)));
    }
    Ok(Design { sample_to_cond, sample_to_covariate, samples_ordered,
                cond1: conditions[0].clone(), cond2: conditions[1].clone() })
}

// ─── Header file parsing (for --no-header) ───────────────────────────────────

/// Read a header file containing column names.
///
/// Accepted formats (auto-detected):
///   • One name per line  (\n-separated)  — used when the file has multiple lines
///   • All names on one line, tab-separated — used when the file is a single
///     non-empty line that contains at least one tab
///
/// In both cases the first name is the feature-ID column; the rest are sample
/// names matched against the design file.
///
/// Returns (feature_col_name, [sample_name, ...]).
fn parse_header_file(path: &PathBuf) -> io::Result<(String, Vec<String>)> {
    let file = File::open(path)
        .map_err(|e| io::Error::new(e.kind(),
            format!("Cannot open header file '{}': {}", path.display(), e)))?;
    let reader = BufReader::new(file);

    // Collect all non-empty lines.
    let lines: Vec<String> = reader.lines()
        .map(|l| l.map(|s| s.trim_end_matches(['\r', '\n']).to_string()))
        .collect::<io::Result<Vec<_>>>()?
        .into_iter()
        .filter(|l| !l.trim().is_empty())
        .collect();

    if lines.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            "Header file is empty."));
    }

    // Decide format:
    //   • Single line with at least one tab  → tab-separated
    //   • Everything else                    → one name per line
    let names: Vec<String> = if lines.len() == 1 && lines[0].contains('\t') {
        // Tab-separated single line
        lines[0].split('\t')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    } else {
        // Newline-separated: each line must be a single name (no internal tabs
        // or spaces, which would indicate a mis-formatted file).
        let mut out = Vec::with_capacity(lines.len());
        for (lineno, line) in lines.iter().enumerate() {
            let name = line.trim();
            if name.contains('\t') || name.contains(' ') {
                return Err(io::Error::new(io::ErrorKind::InvalidData,
                    format!("header file line {}: expected one name per line, \
                             but found whitespace in '{}'", lineno + 1, name)));
            }
            out.push(name.to_string());
        }
        out
    };

    if names.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            "Header file is empty after parsing."));
    }

    let mut names = names;
    let feature_col = names.remove(0);
    if names.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            "Header file must contain at least 2 names \
             (feature-ID column name + at least one sample name)."));
    }
    Ok((feature_col, names))
}

// ─── Normalization file parsing ───────────────────────────────────────────────

fn parse_norm_counts(
    path: &PathBuf,
    required_samples: &[String],
) -> io::Result<HashMap<String, f64>> {
    let file   = File::open(path)?;
    let reader = BufReader::new(file);
    let mut map: HashMap<String, f64> = HashMap::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?.trim().to_string();
        if line.is_empty() { continue; }
        let parts = split_design_line(&line);
        if parts.len() < 2 {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("norm-counts line {}: expected 2 fields, got {}", lineno + 1, parts.len())));
        }
        let sample = parts[0].trim().to_string();
        let sum: f64 = parts[1].trim().parse().map_err(|_|
            io::Error::new(io::ErrorKind::InvalidData,
                format!("norm-counts line {}: cannot parse '{}' as a number", lineno + 1, parts[1])))?;
        if sum <= 0.0 {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("norm-counts line {}: count sum must be > 0, got {}", lineno + 1, sum)));
        }
        map.insert(sample, sum);
    }

    let mut missing: Vec<&str> = required_samples.iter()
        .filter(|s| !map.contains_key(s.as_str()))
        .map(|s| s.as_str()).collect();
    missing.sort();
    if !missing.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("norm-counts file is missing {} sample(s): {:?}", missing.len(), missing)));
    }
    Ok(map)
}

// ─── Per-row processing ───────────────────────────────────────────────────────

/// Everything a worker needs to know — cheaply Arc-shared across threads.
struct WorkCtx {
    idx_cond1:    Vec<usize>,
    idx_cond2:    Vec<usize>,
    pseudo:       f64,
    min_count:    f64,
    max_pvalue:   f64,
    test:         TestType,
    norm_factors: Vec<f64>,
    /// Covariate values aligned to sample column order (same indexing as
    /// idx_cond1/idx_cond2).  Empty when test != Ancova.
    /// covariates[i] is the covariate for the sample at column index i.
    covariates:     Vec<f64>,
    /// When true, store the original input line in HitData for raw output mode.
    keep_raw_line:  bool,
}

struct HitData {
    feature_id: String,
    mean1:  f64,
    mean2:  f64,
    stat1:  f64,   // t / U / t(ANCOVA)
    stat2:  f64,   // df (T-test/ANCOVA) or 0.0 (Wilcoxon)
    raw_p:  f64,
    /// Original input line, stored only when --output-mode raw is active.
    raw_line: Option<String>,
}

impl HitData {
    fn format_line(&self, test: &TestType) -> String {
        let log2fc     = self.mean1 - self.mean2;
        let stat_fields = match test {
            TestType::Ttest    => format!("{:.6}\t{:.2}", self.stat1, self.stat2),
            TestType::Wilcoxon => format!("{:.1}\t.",     self.stat1),
            TestType::Ancova   => format!("{:.6}\t{:.2}", self.stat1, self.stat2),
        };
        format!("{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{:.6e}",
            self.feature_id, self.mean1, self.mean2, log2fc,
            stat_fields, self.raw_p)
    }
}

enum RowResult {
    Hit(HitData),
    FilteredLow,
    FilteredPval,
    NoTest,
}

fn process_row(
    row:  &str,
    ctx:  &WorkCtx,
    buf1: &mut Vec<f64>,
    buf2: &mut Vec<f64>,
) -> RowResult {
    let raw_line = if ctx.keep_raw_line { Some(row.to_string()) } else { None };
    let fields     = split_fields(row);
    let feature_id = fields[0];
    let raw: Vec<f64> = fields[1..].iter()
        .map(|s| s.trim().parse::<f64>().unwrap_or(0.0))
        .collect();

    let get = |i: usize| -> f64 {
        let c = raw.get(i).copied().unwrap_or(0.0);
        if ctx.norm_factors.is_empty() { c }
        else { c / ctx.norm_factors.get(i).copied().unwrap_or(1.0) }
    };

    let all_low = ctx.idx_cond1.iter().chain(ctx.idx_cond2.iter())
        .all(|&i| raw.get(i).copied().unwrap_or(0.0) <= ctx.min_count);
    if all_low { return RowResult::FilteredLow; }

    buf1.clear();
    buf1.extend(ctx.idx_cond1.iter().map(|&i| (get(i) + ctx.pseudo).log2()));
    buf2.clear();
    buf2.extend(ctx.idx_cond2.iter().map(|&i| (get(i) + ctx.pseudo).log2()));

    let mean1 = mean(buf1);
    let mean2 = mean(buf2);

    let result: Option<(f64, f64, f64)> = match ctx.test {
        TestType::Ttest =>
            welch_ttest(buf1, buf2).map(|(t, p, df)| (p, t, df)),
        TestType::Wilcoxon =>
            wilcoxon_ranksum(buf1, buf2).map(|(u, p)| (p, u, 0.0)),
        TestType::Ancova => {
            // Build combined y, group, covariate vectors
            let n1 = ctx.idx_cond1.len();
            let n2 = ctx.idx_cond2.len();
            let n  = n1 + n2;
            let mut y_all   = Vec::with_capacity(n);
            let mut grp_all = Vec::with_capacity(n);
            let mut cov_all = Vec::with_capacity(n);
            for (k, &i) in ctx.idx_cond1.iter().enumerate() {
                y_all.push(buf1[k]);
                grp_all.push(0.0);
                cov_all.push(ctx.covariates.get(i).copied().unwrap_or(0.0));
            }
            for (k, &i) in ctx.idx_cond2.iter().enumerate() {
                y_all.push(buf2[k]);
                grp_all.push(1.0);
                cov_all.push(ctx.covariates.get(i).copied().unwrap_or(0.0));
            }
            ancova(&y_all, &grp_all, &cov_all).map(|(t, p, df)| (p, t, df))
        }
    };

    match result {
        None => RowResult::NoTest,
        Some((p, _, _)) if p > ctx.max_pvalue => RowResult::FilteredPval,
        Some((p, stat1, stat2)) => RowResult::Hit(HitData {
            feature_id: feature_id.to_string(),
            mean1, mean2, stat1, stat2, raw_p: p, raw_line,
        }),
    }
}

// ─── Chunk processing ─────────────────────────────────────────────────────────

struct ChunkOut {
    hits:   Vec<HitData>,
    n_low:  u64,
    n_pval: u64,
    n_na:   u64,
}

fn process_chunk(slab: String, ctx: &WorkCtx) -> ChunkOut {
    let mut out    = Vec::new();
    let mut n_low  = 0u64;
    let mut n_pval = 0u64;
    let mut n_na   = 0u64;
    let mut buf1 = Vec::<f64>::with_capacity(ctx.idx_cond1.len());
    let mut buf2 = Vec::<f64>::with_capacity(ctx.idx_cond2.len());
    for row in slab.lines() {
        if row.is_empty() { continue; }
        match process_row(row, ctx, &mut buf1, &mut buf2) {
            RowResult::Hit(h)       => out.push(h),
            RowResult::FilteredLow  => n_low  += 1,
            RowResult::FilteredPval => n_pval += 1,
            RowResult::NoTest       => n_na   += 1,
        }
    }
    ChunkOut { hits: out, n_low, n_pval, n_na }
}

// ─── Compressed / plain file reader ──────────────────────────────────────────

/// Describes the source of the counts data.
enum CountsSource {
    /// A regular file path (possibly .gz).
    File(PathBuf),
    /// Standard input (already locked / boxed by caller).
    Stdin,
}

fn open_reader_from_source(
    source:  &CountsSource,
    buf_cap: usize,
) -> io::Result<Box<dyn BufRead + Send>> {
    match source {
        CountsSource::Stdin => {
            Ok(Box::new(BufReader::with_capacity(buf_cap, io::stdin())))
        }
        CountsSource::File(path) => {
            open_file_reader(path, buf_cap)
        }
    }
}

fn open_file_reader(path: &PathBuf, buf_cap: usize) -> io::Result<Box<dyn BufRead + Send>> {
    let is_gz = path.extension().map_or(false, |e| e.eq_ignore_ascii_case("gz"));
    if is_gz {
        let child = Command::new("gzip")
            .args(["-dc", &path.to_string_lossy().to_string()])
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .map_err(|e| io::Error::new(e.kind(),
                format!("Failed to spawn zcat for '{}': {} (is zcat installed?)", path.display(), e)))?;
        let stdout = child.stdout.ok_or_else(|| io::Error::new(
            io::ErrorKind::Other, "zcat stdout not captured"))?;
        Ok(Box::new(BufReader::with_capacity(buf_cap, stdout)))
    } else {
        let file = File::open(path)?;
        Ok(Box::new(BufReader::with_capacity(buf_cap, file)))
    }
}

// ─── Pipeline: one streaming pass ────────────────────────────────────────────

fn run_pass<F>(
    pass_label:  &str,
    source:      &CountsSource,
    skip_header: bool,
    ctx:         Arc<WorkCtx>,
    chunk_size:  usize,
    nthreads:    usize,
    mut sink:    F,
) -> io::Result<(u64, u64, u64, u64)>
where
    F: FnMut(&HitData) -> io::Result<()>,
{
    eprintln!("  {}...", pass_label);

    let mut reader = open_reader_from_source(source, chunk_size)?;
    if skip_header {
        let mut header_line = String::new();
        reader.read_line(&mut header_line)?;
    }

    use std::sync::mpsc;
    let (in_tx, in_rx) = mpsc::sync_channel::<(usize, String)>(nthreads * CHANNEL_DEPTH_FACTOR);
    let in_rx = Arc::new(std::sync::Mutex::new(in_rx));
    let (out_tx, out_rx) = mpsc::sync_channel::<(usize, ChunkOut)>(nthreads);

    let reader_thread = thread::spawn(move || {
        let threshold = chunk_size;
        let mut chunk_id = 0usize;
        let mut slab     = String::with_capacity(threshold + 4096);
        let mut line     = String::new();
        loop {
            line.clear();
            match reader.read_line(&mut line) {
                Ok(0) => break,
                Ok(_) => {
                    slab.push_str(&line);
                    if slab.len() >= threshold {
                        if in_tx.send((chunk_id, std::mem::take(&mut slab))).is_err() { break; }
                        chunk_id += 1;
                        slab      = String::with_capacity(threshold + 4096);
                    }
                }
                Err(e) => { eprintln!("Read error: {}", e); break; }
            }
        }
        if !slab.is_empty() { let _ = in_tx.send((chunk_id, slab)); }
    });

    let mut workers = Vec::with_capacity(nthreads);
    for _ in 0..nthreads {
        let irx = Arc::clone(&in_rx);
        let otx = out_tx.clone();
        let c   = Arc::clone(&ctx);
        workers.push(thread::spawn(move || {
            loop {
                let item = irx.lock().unwrap().recv();
                match item {
                    Err(_) => break,
                    Ok((id, slab)) => { let _ = otx.send((id, process_chunk(slab, &c))); }
                }
            }
        }));
    }
    drop(out_tx);
    thread::spawn(move || { reader_thread.join().ok(); });

    let mut pending: BTreeMap<usize, ChunkOut> = BTreeMap::new();
    let mut next   = 0usize;
    let mut n_ok   = 0u64;
    let mut n_low  = 0u64;
    let mut n_pval = 0u64;
    let mut n_na   = 0u64;

    while let Ok((id, chunk_out)) = out_rx.recv() {
        pending.insert(id, chunk_out);
        while let Some(co) = pending.remove(&next) {
            for hit in &co.hits { sink(hit)?; }
            n_ok  += co.hits.len() as u64 + co.n_pval + co.n_na;
            n_low += co.n_low;
            n_pval+= co.n_pval;
            n_na  += co.n_na;
            next  += 1;
        }
    }
    // wait for workers
    for w in workers { w.join().ok(); }
    Ok((n_ok, n_low, n_pval, n_na))
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() -> io::Result<()> {
    let args = Args::parse();

    let nthreads = args.threads
        .unwrap_or_else(|| thread::available_parallelism()
            .map(|n| n.get()).unwrap_or(1))
        .max(1);

    let test_label = match args.test {
        TestType::Ttest    => "Welch T-test",
        TestType::Wilcoxon => "Wilcoxon rank-sum",
        TestType::Ancova   => "ANCOVA (covariate-adjusted T-test)",
    };

    // ── Resolve counts source ──────────────────────────────────────────────
    // --tsv '-'  or  --tsv absent  →  stdin
    // --tsv <path>                 →  file
    let counts_source: CountsSource = match &args.tsv {
        None => CountsSource::Stdin,
        Some(p) if p.to_string_lossy() == "-" => CountsSource::Stdin,
        Some(p) => CountsSource::File(p.clone()),
    };

    // BH requires two passes; stdin can only be read once.
    if args.bh && matches!(counts_source, CountsSource::Stdin) {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
            "--bh (Benjamini-Hochberg) requires two passes and cannot be used \
             when reading from stdin. Please provide a file via --tsv."));
    }

    // 1. Parse design
    eprintln!("Reading design file: {}", args.design.display());
    let design = parse_design(&args.design)
        .map_err(|e| io::Error::new(e.kind(), format!("Design file error: {}", e)))?;
    eprintln!("  Condition 1:    '{}'", design.cond1);
    eprintln!("  Condition 2:    '{}'", design.cond2);
    eprintln!("  Samples:        {}", design.sample_to_cond.len());
    eprintln!("  Test:           {}", test_label);
    eprintln!("  Normalization factor:        {}", args.norm_scale);
    eprintln!("  P-value cutoff: {} ({})",
        args.max_pvalue, if args.bh { "BH q-value" } else { "raw p-value" });
    eprintln!("  BH correction:  {}", if args.bh { "yes" } else { "no" });
    eprintln!("  Threads:        {}", nthreads);
    let chunk_size = IO_BUF_BYTES * args.chunk_factor.max(1);
    eprintln!("  Chunk size:     {} KiB  (factor {})", chunk_size / 1024, args.chunk_factor);

    // Validate covariate presence for ANCOVA
    if matches!(args.test, TestType::Ancova) {
        let missing_cov: Vec<&str> = design.samples_ordered.iter()
            .filter(|s| !design.sample_to_covariate.contains_key(s.as_str()))
            .map(|s| s.as_str())
            .collect();
        if !missing_cov.is_empty() {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("--test ancova requires a numeric third column in the design file, \
                         but {} sample(s) are missing a covariate value: {:?}",
                        missing_cov.len(), missing_cov)));
        }
        eprintln!("  Covariate:      present for all {} samples",
            design.sample_to_covariate.len());
    }

    // 2. Resolve column indices
    //
    //   --no-header <file>  →  read column names from that file; no header
    //                          row is consumed from the counts stream.
    //   (no --no-header)    →  read the first line of the counts stream as
    //                          the header.
    eprintln!("Reading counts from: {}",
        match &counts_source {
            CountsSource::Stdin   => "<stdin>".to_string(),
            CountsSource::File(p) => p.display().to_string(),
        });

    let feature_col:     String;
    let idx_cond1:       Vec<usize>;
    let idx_cond2:       Vec<usize>;
    let has_file_header: bool;           // whether to skip the first data line
    let col_idx_to_name: HashMap<usize, String>;

    match &args.no_header {
        // ── --no-header <HEADER_FILE> ──────────────────────────────────────
        Some(header_path) => {
            eprintln!("  Reading column names from header file: {}",
                header_path.display());
            let (feat_col, sample_names) = parse_header_file(header_path)?;
            feature_col     = feat_col;
            has_file_header = false;   // no header row in counts stream

            // Map sample names to condition using the design file.
            // Column index i corresponds to sample_names[i].
            idx_cond1 = sample_names.iter().enumerate()
                .filter_map(|(i, s)|
                    if design.sample_to_cond.get(s).map_or(false, |c| c == &design.cond1)
                        { Some(i) } else { None })
                .collect();
            idx_cond2 = sample_names.iter().enumerate()
                .filter_map(|(i, s)|
                    if design.sample_to_cond.get(s).map_or(false, |c| c == &design.cond2)
                        { Some(i) } else { None })
                .collect();
            col_idx_to_name = sample_names.iter().enumerate()
                .filter(|(_, s)| design.sample_to_cond.contains_key(s.as_str()))
                .map(|(i, s)| (i, s.clone()))
                .collect();

            let skipped: Vec<&str> = sample_names.iter()
                .filter(|s| !design.sample_to_cond.contains_key(s.as_str()))
                .map(|s| s.as_str())
                .collect();
            eprintln!("  Feature column:   '{}'", feature_col);
            eprintln!("  Samples → '{}': {} (cols {})", design.cond1, idx_cond1.len(),
                idx_cond1.iter().map(|i| (i + 2).to_string()).collect::<Vec<_>>().join(","));
            eprintln!("  Samples → '{}': {} (cols {})", design.cond2, idx_cond2.len(),
                idx_cond2.iter().map(|i| (i + 2).to_string()).collect::<Vec<_>>().join(","));
            if !skipped.is_empty() {
                eprintln!("  Columns in header file but NOT in design (ignored): {:?}", skipped);
            }
        }

        // ── Normal mode: first line of counts stream is the header ─────────
        None => {
            has_file_header = true;

            // We need to peek at the header line before the worker pipeline
            // starts.  For a file we open a fresh reader; for stdin we must
            // read it now (before run_pass opens stdin again — which would
            // give a second BufReader on the same fd, racing with the first).
            //
            // Solution: read the header here via a temporary reader, then
            // signal run_pass to skip the first line (skip_header = true).
            // For stdin this works because the OS fd position advances past
            // the header; run_pass will open a *new* BufReader on stdin and
            // will see only the data lines.

            let header_text = match &counts_source {
                CountsSource::File(path) => {
                    let mut r = open_file_reader(path, chunk_size)?;
                    let mut s = String::new();
                    r.read_line(&mut s)?;
                    s
                }
                CountsSource::Stdin => {
                    // Read header from the real stdin right now.
                    let stdin = io::stdin();
                    let mut r = stdin.lock();
                    let mut s = String::new();
                    r.read_line(&mut s)?;
                    s
                }
            };

            let header_line = header_text.trim_end_matches(['\n', '\r']);
            let headers: Vec<&str> = split_fields(header_line);
            if headers.len() < 2 {
                return Err(io::Error::new(io::ErrorKind::InvalidData,
                    "Header must have at least 2 columns."));
            }
            feature_col = headers[0].to_string();
            let col_cond: Vec<Option<String>> = headers[1..].iter()
                .map(|s| design.sample_to_cond.get(s.trim()).cloned())
                .collect();
            idx_cond1 = col_cond.iter().enumerate()
                .filter_map(|(i, c)| if c.as_deref() == Some(&design.cond1) { Some(i) } else { None })
                .collect();
            idx_cond2 = col_cond.iter().enumerate()
                .filter_map(|(i, c)| if c.as_deref() == Some(&design.cond2) { Some(i) } else { None })
                .collect();
            col_idx_to_name = headers[1..].iter().enumerate()
                .filter(|(_, s)| design.sample_to_cond.contains_key(s.trim()))
                .map(|(i, s)| (i, s.trim().to_string())).collect();
            let skipped: Vec<&str> = headers[1..].iter().zip(col_cond.iter())
                .filter_map(|(s, c)| if c.is_none() { Some(*s) } else { None }).collect();
            eprintln!("  Samples → '{}': {}", design.cond1, idx_cond1.len());
            eprintln!("  Samples → '{}': {}", design.cond2, idx_cond2.len());
            if !skipped.is_empty() {
                eprintln!("  Samples in file but NOT in design (ignored): {:?}", skipped);
            }
        }
    }

    if idx_cond1.is_empty() || idx_cond2.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            "At least one condition has no assigned sample columns."));
    }

    // 3. Build covariate vector aligned to column indices (for ANCOVA)
    let n_sample_cols = idx_cond1.iter().chain(idx_cond2.iter())
        .copied().max().map(|m| m + 1).unwrap_or(0);
    let covariates: Vec<f64> = if matches!(args.test, TestType::Ancova) {
        let mut v = vec![0.0_f64; n_sample_cols];
        for (&idx, name) in &col_idx_to_name {
            if let Some(&cov) = design.sample_to_covariate.get(name.as_str()) {
                v[idx] = cov;
            }
        }
        v
    } else {
        Vec::new()
    };

    // 4. Build normalization factors
    let norm_factors: Vec<f64> = if let Some(ref norm_path) = args.norm_counts {
        eprintln!("Reading norm-counts file: {}", norm_path.display());
        let required: Vec<String> = design.samples_ordered.clone();
        let sum_map = parse_norm_counts(norm_path, &required)
            .map_err(|e| io::Error::new(e.kind(), format!("norm-counts error: {}", e)))?;
        eprintln!("  Loaded sums for {} samples", sum_map.len());
        let mut factors = vec![1.0_f64; n_sample_cols];
        for (&idx, name) in &col_idx_to_name {
            if let Some(&s) = sum_map.get(name.as_str()) {
                factors[idx] = s / args.norm_scale;
            }
        }
        factors
    } else {
        Vec::new()
    };

    // 5. Build shared worker context
    // In BH pass 1, we never need raw lines (we only collect p-values).
    // In single-pass or BH pass 2, keep raw lines only when --output-mode raw.
    let want_raw = matches!(args.output_mode, OutputMode::Raw);

    let ctx_pass1 = Arc::new(WorkCtx {
        idx_cond1:      idx_cond1.clone(),
        idx_cond2:      idx_cond2.clone(),
        pseudo:         args.pseudo,
        min_count:      args.min_count,
        max_pvalue:     if args.bh { 1.0 } else { args.max_pvalue },
        test:           args.test.clone(),
        norm_factors:   norm_factors.clone(),
        covariates:     covariates.clone(),
        keep_raw_line:  want_raw && !args.bh,
    });

    // 6. Open output and write header
    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None    => Box::new(BufWriter::new(stdout.lock())),
    };

    let stat_col = match args.test {
        TestType::Ttest    => "t_stat\tdf",
        TestType::Wilcoxon => "U_stat\t.",
        TestType::Ancova   => "t_stat(ancova)\tdf",
    };
    let adj_col = if args.bh { "\tadj_pvalue" } else { "" };
    // In raw mode we output the original lines with no stats header.
    if !want_raw {
        writeln!(writer,
            "{feature}\tmean_log2_{c1}\tmean_log2_{c2}\tlog2FC({c1}_over_{c2})\t{stat}\tp_value{adj}",
            feature = feature_col, c1 = design.cond1, c2 = design.cond2,
            stat = stat_col, adj = adj_col)?;
    }

    // ── Single-pass ────────────────────────────────────────────────────────
    if !args.bh {
        let (n_ok, n_low, n_pval, n_na) = run_pass(
            "Scanning", &counts_source, has_file_header,
            ctx_pass1, chunk_size, nthreads,
            |hit| {
                if want_raw {
                    if let Some(ref line) = hit.raw_line {
                        writeln!(writer, "{}", line)
                    } else { Ok(()) }
                } else {
                    writeln!(writer, "{}", hit.format_line(&args.test))
                }
            },
        )?;
        eprintln!("Done.");
        eprintln!("  Features processed:              {}", n_ok);
        eprintln!("  Features skipped (low count):    {}", n_low);
        eprintln!("  Features skipped (too few reps): {}", n_na);
        eprintln!("  Features filtered (p > {:.4}):   {}", args.max_pvalue, n_pval);
        eprintln!("  Features reported:               {}", n_ok - n_na - n_pval);
        return Ok(());
    }

    // ── Two-pass BH ────────────────────────────────────────────────────────
    // (only reachable with CountsSource::File, enforced above)
    eprintln!("BH two-pass mode.");
    let mut pval_list: Vec<(String, f64)> = Vec::new();
    let (n_tested, n_low_p1, _, n_na_p1) = run_pass(
        "Pass 1 (computing p-values)", &counts_source, has_file_header,
        Arc::clone(&ctx_pass1), chunk_size, nthreads,
        |hit| { pval_list.push((hit.feature_id.clone(), hit.raw_p)); Ok(()) },
    )?;
    let p1_low = n_low_p1;
    let p1_na  = n_na_p1;
    eprintln!("  Features tested: {}", pval_list.len());

    eprintln!("  Applying Benjamini-Hochberg correction (q ≤ {})...", args.max_pvalue);
    let bh_map: HashMap<String, f64> = benjamini_hochberg(&pval_list, args.max_pvalue);
    eprintln!("  Features passing BH: {}", bh_map.len());

    let passing_ids: HashSet<String> = bh_map.keys().cloned().collect();
    drop(pval_list);

    let ctx_pass2 = Arc::new(WorkCtx {
        idx_cond1:      idx_cond1,
        idx_cond2:      idx_cond2,
        pseudo:         args.pseudo,
        min_count:      args.min_count,
        max_pvalue:     1.0,
        test:           args.test.clone(),
        norm_factors,
        covariates,
        keep_raw_line:  want_raw,
    });

    let passing_ids = Arc::new(passing_ids);
    let bh_map      = Arc::new(bh_map);
    let passing_ids_2 = Arc::clone(&passing_ids);
    let bh_map_2      = Arc::clone(&bh_map);

    run_pass(
        "Pass 2 (writing results)", &counts_source, has_file_header,
        ctx_pass2, chunk_size, nthreads,
        |hit| {
            if passing_ids_2.contains(&hit.feature_id) {
                let adj_p = bh_map_2.get(&hit.feature_id).copied().unwrap_or(1.0);
                if want_raw {
                    if let Some(ref line) = hit.raw_line {
                        writeln!(writer, "{}", line)?;
                    }
                } else {
                    writeln!(writer, "{}\t{:.6e}", hit.format_line(&args.test), adj_p)?;
                }
            }
            Ok(())
        },
    )?;

    eprintln!("Done.");
    eprintln!("  Features processed:              {}", n_tested);
    eprintln!("  Features skipped (low count):    {}", p1_low);
    eprintln!("  Features skipped (too few reps): {}", p1_na);
    eprintln!("  Features reported (BH q ≤ {:.4}): {}", args.max_pvalue, bh_map.len());
    Ok(())
}

// ─── Unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test] fn test_mean_basic() {
        assert!((mean(&[1.0, 2.0, 3.0]) - 2.0).abs() < 1e-10);
    }

    #[test] fn test_sample_variance() {
        let v = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let m = mean(&v);
        assert!((sample_variance(&v, m) - 4.571428571).abs() < 1e-6);
    }

    #[test] fn test_ttest_separated() {
        let a = [10.0_f64, 11.0, 10.5, 10.8, 11.2];
        let b = [1.0_f64,  1.1,  0.9,  1.0,  1.2];
        let (_, p, _) = welch_ttest(&a, &b).unwrap();
        assert!(p < 1e-6, "Expected small p, got {}", p);
    }

    #[test] fn test_ttest_identical() {
        let a = [5.0_f64, 5.0, 5.0];
        let b = [5.0_f64, 5.0, 5.0];
        let (t, p, _) = welch_ttest(&a, &b).unwrap();
        assert!(t.abs() < 1e-10);
        assert!((p - 1.0).abs() < 1e-10);
    }

    #[test] fn test_ttest_insufficient() {
        assert!(welch_ttest(&[1.0], &[2.0, 3.0]).is_none());
    }

    #[test] fn test_ln_gamma_half() {
        let expected = 0.5 * std::f64::consts::PI.ln();
        assert!((ln_gamma(0.5) - expected).abs() < 1e-9);
    }

    #[test] fn test_wilcoxon_separated() {
        let a = [10.0_f64, 11.0, 10.5, 10.8, 11.2];
        let b = [1.0_f64,  1.1,  0.9,  1.0,  1.2];
        let (_, p) = wilcoxon_ranksum(&a, &b).unwrap();
        assert!(p < 0.05, "Expected small p, got {}", p);
    }

    #[test] fn test_wilcoxon_identical() {
        let a = [5.0_f64, 5.0, 5.0, 5.0];
        let b = [5.0_f64, 5.0, 5.0, 5.0];
        let (_, p) = wilcoxon_ranksum(&a, &b).unwrap();
        assert!((p - 1.0).abs() < 1e-6, "Expected p=1, got {}", p);
    }

    #[test] fn test_wilcoxon_insufficient() {
        assert!(wilcoxon_ranksum(&[1.0], &[2.0, 3.0]).is_none());
    }

    #[test] fn test_erfc_known() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-10);
        assert!(erfc(10.0) < 1e-20);
        let p = standard_normal_upper_tail(1.96);
        assert!((p - 0.025).abs() < 0.001, "got {}", p);
    }

    #[test] fn test_bh_basic() {
        let pvals = vec![
            ("geneA".to_string(), 0.001),
            ("geneB".to_string(), 0.008),
            ("geneC".to_string(), 0.039),
            ("geneD".to_string(), 0.041),
            ("geneE".to_string(), 0.850),
        ];
        let result = benjamini_hochberg(&pvals, 0.05);
        assert!(result.contains_key("geneA"), "geneA should pass BH");
        assert!(!result.contains_key("geneE"), "geneE should not pass BH");
        for (gene, raw_p) in &[("geneA", 0.001_f64), ("geneB", 0.008)] {
            if let Some(&adj) = result.get(*gene) {
                assert!(adj >= raw_p - 1e-10,
                    "{}: adj_p {} < raw_p {}", gene, adj, raw_p);
            }
        }
    }

    #[test] fn test_bh_all_pass() {
        let pvals: Vec<(String, f64)> = (1..=5)
            .map(|i| (format!("g{}", i), i as f64 * 0.001)).collect();
        assert_eq!(benjamini_hochberg(&pvals, 0.05).len(), 5);
    }

    #[test] fn test_bh_none_pass() {
        let pvals = vec![("a".to_string(), 0.9), ("b".to_string(), 0.8)];
        assert!(benjamini_hochberg(&pvals, 0.05).is_empty());
    }

    /// ANCOVA test: verify that a covariate explaining all variance gives p≈1,
    /// and that a clear group difference survives covariate adjustment.
    #[test] fn test_ancova_group_effect() {
        // Groups clearly separated, covariate uncorrelated with group
        let y   = [1.0, 1.1, 0.9, 1.0, 5.0, 5.1, 4.9, 5.0_f64];
        let grp = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0_f64];
        let cov = [1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0_f64]; // uncorrelated
        let (_, p, df) = ancova(&y, &grp, &cov).unwrap();
        assert!(p < 0.001, "Expected small p for clear group effect, got {}", p);
        assert!((df - 5.0).abs() < 1e-6, "Expected df=5, got {}", df);
    }

    #[test] fn test_ancova_no_group_effect() {
        // No group effect after covariate adjustment
        let cov = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0_f64];
        // y is purely linear in cov, group adds nothing
        // y is mostly linear in cov but with small noise; group adds nothing
        let y = [3.0, 5.1, 6.9, 9.0, 11.1, 12.9, 15.0, 17.1_f64];
        let grp = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0_f64];
        let (_, p, _) = ancova(&y, &grp, &cov).unwrap();
        assert!(p > 0.5, "Expected large p when group adds nothing, got {}", p);
    }

    #[test] fn test_ancova_insufficient() {
        let y   = [1.0, 2.0, 3.0_f64];
        let grp = [0.0, 0.0, 1.0_f64];
        let cov = [1.0, 2.0, 3.0_f64];
        assert!(ancova(&y, &grp, &cov).is_none());
    }

    #[test] fn test_parse_header_file_newline_separated() {
        use std::io::Write;
        let path = std::env::temp_dir().join("kamscan_test_header_nl.txt");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, "feature_id").unwrap();
            writeln!(f, "sampleA").unwrap();
            writeln!(f, "sampleB").unwrap();
            writeln!(f, "sampleC").unwrap();
        }
        let (feat, samples) = parse_header_file(&path).unwrap();
        assert_eq!(feat, "feature_id");
        assert_eq!(samples, vec!["sampleA", "sampleB", "sampleC"]);
        std::fs::remove_file(&path).ok();
    }

    #[test] fn test_parse_header_file_tab_separated() {
        use std::io::Write;
        let path = std::env::temp_dir().join("kamscan_test_header_tab.txt");
        {
            let mut f = File::create(&path).unwrap();
            write!(f, "feature_id\tsampleA\tsampleB\tsampleC").unwrap();
        }
        let (feat, samples) = parse_header_file(&path).unwrap();
        assert_eq!(feat, "feature_id");
        assert_eq!(samples, vec!["sampleA", "sampleB", "sampleC"]);
        std::fs::remove_file(&path).ok();
    }
}
