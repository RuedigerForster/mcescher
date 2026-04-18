# Discussion

## Open issues

- peak detected at the end of inhibit period 1
- areas peak widths look strange in peak detection
- uncompress baseline and chromatogram before integration
- yaml configuration file
  - inhibit begin end
  - minimum peak width
  - peak names RT

---

## QuanFormer — evaluation (2026-04-11)

Repository: `./QuanFormer`

### What it is

QuanFormer is a Python/PyTorch all-in-one pipeline for smoothing, baseline
estimation, peak detection and integration of LC-MS chromatograms.  The central
idea is to convert each Extracted Ion Chromatogram (EIC) to a JPEG image and
run a DETR object detector (ResNet50 backbone + Transformer encoder/decoder)
that predicts bounding boxes around peaks.  Peak areas are then computed by
trapezoidal integration within the detected boundaries.

### Architecture

```
mzML  →  EIC extraction (PPM-tolerance)  →  JPEG rendering
      →  DETR inference (checkpoint0029.pth)
      →  bbox → RT coordinate conversion  →  trapz integration
      →  post-processing (duplicate removal, pivot to compound × sample)
```

Two modes: targeted (compound list with expected m/z and RT) and untargeted
(XCMS centWave via an R subprocess).

### Strengths

- **DETR detector is genuinely good.** Trained on real chromatography data; handles
  overlapping, asymmetric and co-eluting peaks better than derivative-based
  threshold detectors.  Confidence scores provide a natural quality gate.
- Clean separation of concerns across utility modules.
- GUI uses QThread correctly — non-blocking for long operations.
- DETR model is lightweight (1 encoder + 1 decoder layer); inference is fast
  even on CPU.
- Obiwarp RT alignment included via XCMS for untargeted data.

### Critical weaknesses

1. **Hardcoded image geometry in quantification** (`quantify.py`).
   The pixel-to-RT coordinate inversion relies on magic constants (50, 400, 30,
   320) that encode the matplotlib figure padding.  Any change to the
   visualisation code silently produces wrong areas.  This is the most
   fundamental flaw for a quantification tool.

2. **Maximum 3 peaks per ROI.**  `num_queries=3` in the DETR config.  Additional
   peaks are silently dropped with no warning.

3. **R subprocess via `os.system()`** — no return code checked, no path escaping;
   spaces in filenames break it silently.

4. **No detection observability.**  Failed detections (empty scores) return zero
   area; indistinguishable from a true absence of the compound.

5. **Pickle-based pipeline state** — fragile across code changes; no streaming.

6. **Old, pinned dependencies** (PyTorch 1.13.1, numpy 1.24.4).

### Fit with the current pipeline

QuanFormer and this pipeline address the same problem from opposite directions:

| Aspect            | This pipeline                        | QuanFormer              |
|-------------------|--------------------------------------|-------------------------|
| Baseline          | Algorithmic ensemble (arPLS, BEADS…) | Implicit in model       |
| Peak detection    | Derivative zero-crossing             | DETR on EIC image       |
| Integration       | PD / TS / Gauss / EGH                | Trapezoid               |
| Uncertainty       | CV% across ensemble configs          | None                    |
| Explainability    | Every step auditable                 | Black box               |

QuanFormer is not suitable as a drop-in replacement.  It provides no uncertainty
estimate and no interpretable baseline, which are both requirements for
regulatory-grade quantification.

### Parts worth extracting

**High value — DETR as a smarter `find_peaks()` front-end.**
The trained detector is the strongest part of the repo.  A viable integration
path:

1. Write signal + RT to a temp CSV from R.
2. Python script renders each signal as a JPEG with *controlled, recorded*
   axis limits; runs DETR; maps the normalised bounding box output
   `[cx, cy, w, h] ∈ [0,1]` directly to RT boundaries via:

   ```
   RT_left  = xlim[0] + (cx - w/2) * (xlim[1] - xlim[0])
   RT_right = xlim[0] + (cx + w/2) * (xlim[1] - xlim[0])
   ```

   This completely avoids QuanFormer's pixel-geometry fragility.
3. R wrapper converts bbox output to `find_peaks()` format
   (`idx`, `RT`, `height`, `width`, `area`, `lo_idx`, `hi_idx`).
4. `pipeline_summary()` gets a `peak_finder` parameter (`"derivative"` or
   `"detr"`); everything downstream (ensemble, integration, LOD/LOQ) is
   unchanged.

**Medium value — EIC extraction (`extract_eic.py`).**
Clean PPM-toleranced EIC builder with joblib parallelism.  Worth adapting if
the pipeline moves to direct mzML ingestion.

**Low value — XCMS wrapper pattern.**
The R-subprocess approach (`detect_helper.py` + `find_peaks.R`) is a useful
pattern but needs hardening (`subprocess.run(..., check=True)`, proper path
quoting).

### Status / revised placement

**Not to be integrated here.**  QuanFormer fits better into the second pipeline
approach: CNN-based learning of the analyst's integration style.

In that context the image-based representation is an asset rather than a
workaround — a CNN processes the EIC image the same way a trained analyst does,
which is precisely the goal when the objective is to replicate human integration
decisions.  The DETR bounding boxes become training labels or a prior rather
than a final quantification result, making the pixel-geometry fragility in
`quantify.py` irrelevant.

The DETR detector (`checkpoint0029.pth`) and the EIC rendering pipeline
(`extract_eic.py`, `plot_utils.py`) are the most directly reusable components
for that work.
