# m.c.escher

<p align="center">
  <img src="logo.png" width="160" alt="m.c.escher logo"/>
</p>

> Automated chromatographic peak integration with uncertainty quantification.

---

## Motivation

Manual peak integration is the last routine step in chromatographic data analysis
that still depends heavily on expert judgement. Results vary between operators,
software versions, and integration parameter sets — yet the reported area is
usually presented as a single number with no uncertainty attached.

**m.c.escher** addresses this by treating baseline placement as an ensemble
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
- **Multi-CDS import** via [chromConverter](https://cran.r-project.org/package=chromConverter)
  (CDF/AIA, Waters ARW, Agilent, Shimadzu, and others)
- **LOD / LOQ** — ICH Q2(R1) signal-to-noise estimation per chromatogram

---

## Status

> Early development. API is not yet stable.

The pipeline runs end-to-end on single-channel GC-FID/TCD data in CDF/AIA format.
The ML integration model (the core of the package name) is in the design phase;
the current release focuses on the classical baseline-ensemble pipeline that will
generate training labels.

---

## Installation

```r
# From GitHub (development version)
# install.packages("remotes")
remotes::install_github("RuedigerForster/Integration")
```

CRAN submission is planned once the API stabilises.

---

## Quick start

```r
library(mcescher)   # package name on CRAN will be mcescher

# Import CDF files from a directory
chromas <- import_chromas("path/to/data/", pattern = "[.]cdf$")

# Run the full pipeline on one chromatogram
result <- pipeline_summary(
  signal_mat  = chromas,
  method_yaml = "method.yaml"
)

# Peak table with areas and expanded uncertainties
print(result$peaks)
```

---

## Pipeline overview

```
CDF / ARW file
      │
      ▼
 chromConverter::read_cdf()          # multi-CDS import
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

The pipeline is controlled by a YAML file:

```yaml
method:
  name: "GC-FID light hydrocarbons"
  RT_unit: min

inhibit:
  - [0.00, 1.10]    # solvent front
  - [20.38, 27.00]  # tail region

peaks:
  - name: Methane
    RT_ref:    1.37
    RT_window: 0.20

processing:
  x_factor: 50      # compression factor
  denoise:   true

integration:
  methods: [PD, TS, gauss, EGH]
  run_ensemble: true
  ensemble_method: TS
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
  title   = {m.c.escher: Automated Chromatographic Peak Integration
             with Uncertainty Quantification},
  year    = {2026},
  url     = {https://github.com/RuedigerForster/Integration}
}
```
