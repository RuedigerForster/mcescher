# mcescher

<p align="center">
  <img src="man/figures/mcescher.png" width="160" alt="mcescher logo"/>
</p>

> Automated chromatographic peak integration with GUM-compliant uncertainty quantification.

---

## For ML researchers

Gas chromatography produces a signal trace over time. Compounds elute as peaks;
their areas encode concentration. Training a neural network to integrate these
peaks requires reliable (signal, area) pairs — but most GC software ships areas
with no indication of whether the integration was clean or marginal.

**mcescher** runs 15 baseline algorithms in parallel, selects a consensus by
mutual agreement, and attaches a three-component GUM uncertainty to every area.
That uncertainty is a direct, quantitative label-quality signal: a 0.2 % CV is
a clean label; a 40 % CV is noise.

### Get a training set in three steps

```r
# install.packages(c("remotes", "nanoparquet", "chromConverter"))
# remotes::install_github("RuedigerForster/mcescher")

library(mcescher)

# 1. Run mcescher on your data (see §Running mcescher below)
#    peaks_raw <- peaks(pipeline_summary(...))

# 2. Apply quality gates
peaks_ok <- peaks_raw |>
  replicate_quality_gate() |>   # remove within-sequence outliers
  loq_filter()             |>   # remove below-LOQ labels
  (\(df) df[df$quality_ok, ])()

# 3. Export paired (signal window, area, uncertainty) to Parquet
export_training_set(
  peaks_df       = peaks_ok,
  data_root      = "path/to/chromatogram/files",
  output_parquet = "training_set.parquet"
)
```

Read in Python:

```python
import pandas as pd, numpy as np

df = pd.read_parquet("training_set.parquet")
X  = np.stack(df["signal"])       # (N, 512) float32 — normalised signal window
y  = df["area_cons"].values        # (N,)  area in mV·s
u  = df["u_cons"].values           # (N,)  GUM k=1 standard uncertainty
cv = df["area_cv_pct"].values      # (N,)  relative uncertainty in %
```

### Understanding the quality columns

| Column | What it means for training |
|---|---|
| `area_cv_pct` | Relative standard uncertainty of the label (%). Filter on this. Labels with CV > 5 % are near the detection limit — the area estimate is dominated by baseline noise, not peak signal. |
| `u_cons` | Absolute standard uncertainty (mV·s, k=1). Use as per-sample loss weight: `loss = ((ŷ - y) / u_cons)²` gives higher weight to clean labels. |
| `Z_spread` | Peak difficulty score. High values indicate fused peaks, sloping baseline, or other integration challenges. |
| `aic_winner` | Integration method selected by AICc (PD / TS / gauss / EGH). Consistent method across replicates is a sign of a well-resolved peak. |

### Quality gate parameters

```r
replicate_quality_gate(
  peaks_df,
  min_n          = 5,    # minimum replicates before gate activates
  z_threshold    = 3.5,  # Iglewicz-Hoaglin modified Z-score cutoff
  max_seq_cv_pct = 3.0   # flag entire compound if sequence CV% exceeds this
)

loq_filter(
  peaks_df,
  max_area_cv_pct = 5.0  # drop rows where area_cv_pct exceeds this
)
```

---

## For chromatographers

Manual peak integration still depends heavily on expert judgement. Results vary
between operators, software versions, and integration parameter sets — yet the
reported area is almost always presented as a single number with no uncertainty
attached. Accreditation bodies accept this. Downstream statistics ignore it.
The area walks into calibration as if it were exact.

**mcescher** treats baseline placement as what it actually is: an estimation
problem with multiple plausible solutions. Fifteen baseline configurations run
in parallel; their results are combined by correlation-weighted consensus;
and a three-component GUM uncertainty is attached to every area — baseline
ambiguity, detector noise propagated through integration, and denoising
distortion, each reported separately and combined in quadrature.

The second step is the one most integrators miss. Within a batch, all peaks
share the same injection — so their areas are correlated. A baseline that
wanders under one peak is likely wandering under adjacent peaks too. mcescher
exploits this: Z-scores are computed relative to the expected smooth area curve
across the chromatogram, and peaks with deviating Z-scores have their
uncertainty inflated accordingly. The `Z_spread` column is the integration
difficulty indicator that CDS software has never reported.

The result is an area estimate with a realistic, traceable uncertainty — not
because it makes the number look more scientific, but because it tells you
which areas to trust and which to treat with caution.

---

## Reference dataset

The OGE-UDE Grob mix dataset (2693 GC-FID replicate injections, Agilent
ChemStation .D format) used to validate mcescher and generate training labels
for the **stanislaw** CNN integrator is published on Zenodo:

> Görs, P. E., Schmitz, O., Forster, R. (2026).
> *Grob Mix GC-FID Replicate Injections for Atmospheric Pressure Peak Area
> Correction (APPAC) and Integration Benchmarking.*
> Zenodo. https://doi.org/10.5281/zenodo.19946728

This dataset is not bundled with the package. Download and unpack it, then
point `data_root` in `export_training_set()` at the extracted directory.

---

## Installation

```r
remotes::install_github("RuedigerForster/mcescher")
```

Optional dependencies (install as needed):

```r
install.packages(c(
  "nanoparquet",    # Parquet export via export_training_set()
  "chromConverter", # reading Agilent .D directories
  "ncdf4"           # reading AIA/CDF files
))
```

CRAN submission planned once the API stabilises.

---

## Running mcescher

The pipeline takes a numeric matrix as its primary input.

```r
library(mcescher)

cfg    <- read_method_config("method.yaml")
result <- pipeline_summary(
  signal_mat   = sig,      # n_samples × n_chroms numeric matrix
  baseline_mat = bl,       # arPLS baseline matrix (same dimensions)
  RT           = RT,       # retention time vector (minutes)
  method_config = cfg,
  run_ensemble  = TRUE
)

peaks(result)              # data frame of areas, uncertainties, LOD/LOQ
plot(result, chrom = 1)    # signal + consensus baseline + ±2σ band
```

For full-batch processing of Agilent .D directories see
`R/orphans/run_ude.R`, which handles ATSA alignment, SASS denoising, baseline
estimation, integration, and quality gating for an entire dataset in one call.

---

## How it works

```
n_samples × n_chroms signal matrix
      │
      ▼
 Block-average compression (×50)
      │
      ▼
 ATSA alignment                  — correct retention-time drift across runs
      │
      ▼
 SASS denoising                  — sparse signal smoothing (L1)
      │  └─ save pre-SASS ──────────────────────────────────► u_sass
      ▼
 arPLS baseline
      │
      ▼
 Decompress → full resolution
      │
      ▼
 Peak detection + integration    — PD / TS / Gaussian / EGH
      │
      ▼
 Baseline ensemble (15 configs)  — airPLS · arPLS · LMV · FlatFit · BEADS
      │  ├─ Step 1: correlation weights ───────────────────► area_cons
      │  └─ Step 2: Z-score penalty ───────────────────────► u_baseline
      │
      ▼
 Noise estimation ───────────────────────────────────────► u_noise
      │
      ▼
 u_area = √(u_baseline² + u_noise² + u_sass²)
      │
      ▼
 replicate_quality_gate() → loq_filter() → export_training_set()
```

---

## Uncertainty model

| Component | Symbol | Source |
|---|---|---|
| Baseline placement ambiguity | `u_cons` | Ensemble spread + Z-score inflation |
| Detector noise propagated through integration | `u_noise` | Noise floor × √(peak width × Δt) |
| SASS denoising distortion | `u_sass` | \|∫(pre − post SASS)\| over peak window |

Combined: **U = 2 · √(u_baseline² + u_noise² + u_sass²)** (k = 2, ~95 %).

---

## Output columns

| Column | Description |
|---|---|
| `area_cons` | Consensus area — best point estimate (mV·s) |
| `u_cons` | Baseline standard uncertainty (k=1) |
| `u_noise` | Noise-propagation standard uncertainty |
| `u_sass` | SASS-distortion standard uncertainty |
| `u_area` | Combined standard uncertainty (k=1) |
| `U_area` | Expanded uncertainty (k=2, ~95 %) |
| `area_cv_pct` | Relative standard uncertainty (%) — primary label-quality indicator |
| `Z_spread` | Weighted mean \|Z\| — integration difficulty |
| `penalty` | Step-2 inflation factor (1 = clean, > 1 = ambiguous baseline) |
| `aic_winner` | AICc-selected integration method |
| `area_PD` | Area by Perpendicular Drop |
| `area_TS` | Area by Tangent Skim |
| `area_gauss` | Area by Gaussian fit |
| `area_EGH` | Area by EGH fit |

---

## Supported baseline algorithms

| Algorithm | Reference | License |
|---|---|---|
| airPLS | Zhang et al., *Analyst* 2010 | BSD-3 |
| arPLS | Baek et al., *Analyst* 2015 | BSD-3 |
| LMV-RSA | — | BSD-3 |
| FlatFit | MOCCA, Bayer AG | MIT |
| BEADS | Ning, Selesnick & Duval, *Chemom. Intell. Lab. Syst.* 2014 | — |

---

## Name and logo

The package name honours **M. C. Escher** (1898–1972), who was a devoted admirer
of the Amalfi coast. The landscape — its staggered cliffs, terraced villages,
and interlocking planes of rock and sea — left a lasting mark on his visual
language. The idea of decomposing a chromatogram into interlocking,
self-consistent baseline segments traces the same geometric intuition.

The logo is a photograph of those same cliffs near Maiori, rendered in shades
of blue, taken at the moment the project was conceived.

---

## License

GPL (≥ 3)

---

## Citation

```bibtex
@software{forster2026mcescher,
  author  = {Forster, Rüdiger},
  title   = {mcescher: Automated Chromatographic Peak Integration
             with Uncertainty Quantification for Target Compound Analysis},
  year    = {2026},
  url     = {https://github.com/RuedigerForster/mcescher}
}
```
