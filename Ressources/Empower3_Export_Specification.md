# Empower 3 Export Specification — ML Training Data

## Overview
Two Export Methods are needed per batch, applied to the same set of results:
1. **Raw chromatogram data** → CDF files (one per channel per injection)
2. **Peak table** → ASCII file (one per result)

Both are linked by the Result ID that Empower appends automatically to filenames.

---

## Export Method 1 — Raw Chromatogram Data

**Raw Data tab:**
| Setting | Value |
|---|---|
| Format | AIA |
| Destination | File |
| Root name | `raw_` |

Produces: `raw_<ChannelID>_<ResultID>.cdf`

---

## Export Method 2 — Peak Table (ASCII)

**Fields tab:**
| Setting | Value |
|---|---|
| Export ASCII File | Checked |
| Export Table Data | Checked |
| Export Field Data | Checked |
| Delimiter | Tab or comma (consistent) |
| Root name | `peaks_` |

**Fields to include — confirm exact names in your Empower system:**

| Field | Purpose |
|---|---|
| Sample Name | Links injection to context |
| Vial | Identifies sample position |
| Injection Number | Replicate identifier |
| Channel Name | Which of the 2 channels |
| Peak Name | Compound identity |
| Peak Number | Order in chromatogram |
| Retention Time | Peak apex |
| Start | Peak start time → label |
| End | Peak end time → label |
| Area | Integration result |
| Height | Peak height |
| Baseline Code | Integration type → label |
| Manual Integration | Expert correction flag → label |
| Signal to Noise | Quality filter |
| Width at Half Height | Window sizing |

Produces: `peaks_<ResultID>.txt`

---

## Storage
- Estimate ~20–50 GB for a full batch of 5000 runs
- Organise by batch in subfolders: `Batch01/`, `Batch02/`, etc.
- Keep raw CDF and peak tables in separate subfolders per batch:
  ```
  Batch01/
    raw/      <- .cdf files
    peaks/    <- .txt files
  ```

---

## Batch Export Procedure
1. Open the Results tab in Empower
2. Filter to the target batch
3. Select all results (Ctrl+A)
4. Right-click → Export → select Export Method 1 → OK
5. Wait for background processing to complete
6. Repeat with Export Method 2
7. Verify file counts: should be `n_injections × 2` CDF files and `n_injections` peak tables

---

## Pilot First
Before exporting a full batch, run both Export Methods on **10–20 results** and verify in R that `chromConverter` reads the CDF files and that the peak table fields are complete and correctly named. Adjust field selection before committing to the full export.
