# mcescher

<p align="center">
  <img src="man/figures/mcescher.png" width="160" alt="mcescher logo"/>
</p>

> Automated chromatographic peak integration with uncertainty quantification for target compound analysis.

---

## Motivation

Manual peak integration is a crucial step in chromatographic data analysis
that still depends heavily on expert judgement. Results vary between operators,
software versions, and integration parameter sets — yet the reported area is
usually presented as a single number with no uncertainty attached.

**mcescher** addresses this by treating baseline placement as an ensemble
estimation problem: multiple baseline algorithms and parameter combinations are
run in parallel, their results are combined using a correlation-weighted consensus,
and a realistic uncertainty is derived from three independent sources:

| Component | Symbol | Source |
|---|---|---|
| Baseline placement ambiguity | u_baseline | Two-step ensemble (correlation weights + Z-score penalty) |
| Detector noise propagated through integration | u_noise | Noise floor × √(peak width × Δt) |
| SASS denoising distortion | u_sass | \|∫(pre − post SASS)\| over peak window |

These combine in quadrature: **U = 2 · √(u_baseline² + u_noise² + u_sass²)** (k = 2, ~95 %).

---

## Premises

Current instrument technology is highly advanced, so that the sampled signal trace can be considered artifact-free in most cases. The dominant uncertainty sources of the peak areas are the injected sample amounts — which, because all peaks share the same injection, introduce a high degree of correlation across the chromatogram — and the subjective estimation of the baseline course. Consequently, a baseline that preserves this inter-peak correlation is preferred as the consensus estimate.

We expect the ideal peak areas to be affected only by Gaussian noise, and their uncertainties to follow a smooth curve (dominated by detector non-linearity). Occasionally, baseline drift and wandering may disturb individual peaks and break this pattern. Therefore, in a second step, Z-scores of the peak areas are calculated relative to the expected smooth curve, and peaks with deviating Z-scores have their uncertainty inflated accordingly.

---

## Features

- **Multi-algorithm baseline ensemble** — airPLS, arPLS, LMV-RSA, FlatFit, BEADS;
  15 parameter combinations per run
- **Two-step consensus**
  - *Step 1*: correlation-based config weights — globally rogue baselines suppressed
  - *Step 2*: Z-score penalty — locally difficult peaks (fused, sloping baseline) flagged
    and their uncertainty inflated
- **Three-component uncertainty** — baseline ambiguity, noise, and denoising distortion
  reported separately and combined
- **Peak integration methods** — Perpendicular Drop, Tangent Skim, Gaussian fit, EGH fit
- **Time-shift alignment** — ATSA (automatic time-shift alignment) across sample batches
- **LOD / LOQ** — ICH Q2(R1) signal-to-noise estimation per chromatogram

---

## Status

> Early development. API is not yet stable.

The package accepts a pre-imported numeric matrix as input (see Quick start below).
Chromatogram import from CDF/AIA or vendor-specific formats is not bundled —
use [chromConverter](https://cran.r-project.org/package=chromConverter) or
similar for file reading, then pass the resulting matrix to `pipeline_summary()`.

The ML integration model (the core of the package name) is in the design phase;
the current release focuses on the classical baseline-ensemble pipeline that will
generate training labels.

---

## Installation

```r
# From GitHub (development version)
# install.packages("remotes")
remotes::install_github("RuedigerForster/mcescher")
```

CRAN submission is planned once the API stabilises.

---

## Quick start

The pipeline takes a numeric matrix as its primary input:

| Argument | Type | Description |
|---|---|---|
| `signal_mat` | `n_samples × n_chroms` matrix | Aligned, SASS-denoised signal |
| `baseline_mat` | `n_samples × n_chroms` matrix | Primary arPLS baseline (NA in inhibited regions) |
| `RT` | numeric vector (length `n_samples`) | Retention time in minutes |

`signal_mat` must have a `retention_time` attribute (minutes) set on the column dimension, or `RT` must be passed explicitly.

```r
library(mcescher)

# Assume `sig`, `bl`, and `RT` are prepared from your import step
cfg    <- read_method_config("method.yaml")
result <- pipeline_summary(
  signal_mat   = sig,
  baseline_mat = bl,
  RT           = RT,
  method_config = cfg,
  run_ensemble  = TRUE
)

# ChromResult: peaks with areas, uncertainties, LOD/LOQ
print(result)
as.data.frame(result)      # the peaks table
plot(result, chrom = 1)    # signal + consensus baseline + ±2σ band
```

---

## Pipeline overview

```
n_samples × n_chroms signal matrix  (+ RT vector)
      │
      ▼
 Block-average compression (×50)     # reduce to ~20 Hz working resolution
      │
      ▼
 ATSA alignment                      # correct retention-time drift across runs
      │
      ▼
 SASS denoising                      # sparse signal denoising (L1 regularisation)
      │  └─ save pre-SASS signal ────────────────────────► u_sass
      ▼
 arPLS baseline                      # primary baseline for peak detection
      │
      ▼
 Decompress to full resolution
      │
      ▼
 Peak detection                      # smoothed derivative + amplitude threshold
      │
      ▼
 Peak integration                    # PD / TS / Gaussian / EGH
      │
      ▼
 Baseline ensemble (15 configs)      # airPLS · arPLS · LMV · FlatFit · BEADS
      │  ├─ Step 1: correlation weights ──────────────────► area_cons
      │  └─ Step 2: Z-score penalty ──────────────────────► u_baseline
      │
      ▼
 Noise estimation (LOD / LOQ) ──────────────────────────► u_noise
      │
      ▼
 u_area = √(u_baseline² + u_noise² + u_sass²)
      │
      ▼
 ChromResult (peaks · uncertainties · LOD/LOQ · alignment)
```

---

## Output columns

| Column | Description |
|---|---|
| `area_cons` | Correlation-weighted consensus area (best point estimate) |
| `u_step1` | Baseline uncertainty before Z-score penalty |
| `Z_spread` | Weighted mean \|Z\| — integration difficulty indicator |
| `penalty` | Step-2 inflation factor (1 = clean peak, > 1 = ambiguous) |
| `u_cons` | Penalised baseline standard uncertainty |
| `u_noise` | Noise-propagation standard uncertainty |
| `u_sass` | SASS-distortion standard uncertainty |
| `u_area` | Combined standard uncertainty (k = 1) |
| `U_area` | Expanded uncertainty (k = 2, ~95 %) |
| `area_PD` | Area by Perpendicular Drop |
| `area_TS` | Area by Tangent Skim |
| `area_gauss` | Area by Gaussian fit |
| `area_EGH` | Area by Exponentially Modified Gaussian fit |

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

## Method configuration

The pipeline is controlled by a YAML file passed to `read_method_config()`:

```yaml
method:
  name: "GC light hydrocarbons"
  RT_unit: min          # "min" scales areas to counts·s (×60)

inhibit:
  - [0.00, 1.10]        # valve-switching artifacts, solvent front
  - [20.38, 27.00]      # tail region

# Baseline-inhibit: RT windows replaced by linear interpolation in the baseline.
# Useful when a large spike would distort the ensemble algorithms.
baseline_inhibit:
  - [0.00, 1.20]

peaks:
  - name: Methane
    RT_ref:    1.37
    RT_window: 0.20
    # merge_window: merge all detected sub-peaks within ±N min of RT_ref into one.
    # merge_window: 0.15
    # ensemble_method: override the ensemble consensus with a single-method area.
    # ensemble_method: TS

processing:
  x_factor: 50          # block-average compression factor
  denoise:   true

# Global denoising parameters (SASS):
denoising:
  d:   1
  fc:  0.011
  K:   1
  lam: 0.2

# Per-segment denoising overrides (e.g. solvent front needs different fc):
denoising_segments:
  - rt_min: 0.0
    rt_max: 2.0
    fc:  0.03
    lam: 0.5

baseline:
  lambda: 1.0e7
  ratio:  1.0e-6

detection:
  amp_thresh:          0
  smooth_width_factor: 3

integration:
  methods: [PD, TS, gauss, EGH]
  run_ensemble:    true
  ensemble_method: TS

# debug_baselines: save a PNG of all ensemble baselines for this peak.
# peaks:
#   - name: PropanePeak
#     debug_baselines: true
```

---

## Name and logo

The package name honours **M. C. Escher** (1898–1972), who was a devoted admirer
of the Amalfi coast and visited it repeatedly between 1925 and 1936. The
landscape — its staggered cliffs, terraced villages, and interlocking planes of
rock and sea — left a lasting mark on his visual language.

The logo is a photograph of those same cliffs near Maiori, rendered in shades of
blue, taken at the moment the project was conceived. The connection is not merely
wordplay: the idea of decomposing a chromatogram into interlocking, self-consistent
segments traces the same geometric intuition that runs through Escher's work.

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
