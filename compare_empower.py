"""
compare_empower.py — compare mcescher integration vs Empower (CDS) peak table.

Run after compare_empower.R has written mcescher_peaks.csv.

Usage:
    python compare_empower.py <data_dir>
    e.g. python compare_empower.py data/UNIS_FID_9000002_908779
"""

import sys
import os
import re
import numpy as np
import pandas as pd
import netCDF4 as nc

# ── Empower peak table extraction ─────────────────────────────────────────────

def extract_empower_peaks(cdf_path: str) -> pd.DataFrame:
    """Read Empower CDS peak table from a CDF/AIA file."""
    ds = nc.Dataset(cdf_path)
    rt     = np.array(ds.variables["peak_retention_time"][:])
    area   = np.array(ds.variables["peak_area"][:])
    height = np.array(ds.variables["peak_height"][:])
    amount = np.array(ds.variables["peak_amount"][:])
    raw_names = ds.variables["peak_name"][:]
    names = []
    for row in raw_names:
        s = "".join(b.decode("latin-1") for b in row if b != b"\x00")
        # strip trailing non-ASCII Empower metadata bytes
        s = re.sub(r"[^\x20-\x7E].*", "", s).strip()
        names.append(s)
    ds.close()

    SENTINEL = 1e30
    area   = np.where(area   > SENTINEL, np.nan, area)
    height = np.where(height > SENTINEL, np.nan, height)
    amount = np.where(amount > SENTINEL, np.nan, amount)

    mask = rt > 0
    df = pd.DataFrame({
        "empower_name":   np.array(names)[mask],
        "empower_RT_s":   rt[mask],
        "empower_RT_min": rt[mask] / 60.0,
        "empower_area":   area[mask],
        "empower_height": height[mask],
        "empower_amount": amount[mask],
    })
    df["source_file"] = os.path.basename(cdf_path)
    return df


def load_empower_all(data_dir: str) -> pd.DataFrame:
    files = sorted(f for f in os.listdir(data_dir) if f.lower().endswith(".cdf"))
    frames = []
    for f in files:
        try:
            df = extract_empower_peaks(os.path.join(data_dir, f))
            frames.append(df)
        except Exception as e:
            print(f"  Warning: could not read {f}: {e}")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# ── Comparison ────────────────────────────────────────────────────────────────

NAMED_PEAKS = ["Methane", "Ethane", "Propane", "iso-Butane", "n-Butane"]
EMPOWER_MAP = {
    "CH4":       "Methane",
    "C2H6":      "Ethane",
    "C3H8":      "Propane",
    "i-C4H10":   "iso-Butane",
    "i-C4H10z":  "iso-Butane",
    "n-C4H10":   "n-Butane",
    "n-C4H10z":  "n-Butane",
}


def main(data_dir: str):
    mcescher_csv = os.path.join(data_dir, "mcescher_peaks.csv")
    if not os.path.exists(mcescher_csv):
        print(f"ERROR: {mcescher_csv} not found — run compare_empower.R first")
        sys.exit(1)

    mc  = pd.read_csv(mcescher_csv)
    emp = load_empower_all(data_dir)

    # Normalise Empower names using known mapping
    emp["peak_name"] = emp["empower_name"].map(EMPOWER_MAP)

    # Keep only the five named peaks
    mc_named  = mc[mc["name"].isin(NAMED_PEAKS)].copy()
    emp_named = emp[emp["peak_name"].isin(NAMED_PEAKS)].copy()

    # Merge on chromatogram file + peak name
    # mcescher uses full source file name in 'chrom'; Empower uses basename
    mc_named["file_base"] = mc_named["chrom"].apply(os.path.basename)
    emp_named["file_base"] = emp_named["source_file"]

    merged = pd.merge(
        mc_named[["file_base", "name", "RT", "area_cons", "U_area",
                   "area_PD", "area_TS", "area_gauss", "area_EGH",
                   "area_AIC", "aic_winner",
                   "u_noise", "u_sass", "area_cv_pct"]],
        emp_named[["file_base", "peak_name", "empower_RT_min",
                   "empower_area", "empower_height", "empower_amount"]],
        left_on=["file_base", "name"],
        right_on=["file_base", "peak_name"],
        how="outer",
    )

    merged["area_ratio"] = merged["area_cons"] / merged["empower_area"]
    merged["area_diff_pct"] = (merged["area_cons"] - merged["empower_area"]) \
                               / merged["empower_area"] * 100
    merged["aic_ratio"]    = merged["area_AIC"] / merged["empower_area"]
    merged["aic_diff_pct"] = (merged["area_AIC"] - merged["empower_area"]) \
                               / merged["empower_area"] * 100

    # ── Summary table per peak — ensemble (area_cons) ─────────────────────────
    print("\n═══ mcescher vs Empower — area_cons (baseline ensemble) ═══\n")
    print(f"{'Peak':<12} {'n':>4}  {'Empower mean':>14}  "
          f"{'mcescher mean':>14}  {'Ratio':>7}  {'Diff%':>8}  {'U_area%':>8}  {'CV%':>7}")
    print("─" * 85)

    for peak in NAMED_PEAKS:
        sub = merged[merged["name"] == peak].dropna(subset=["area_cons", "empower_area"])
        if sub.empty:
            print(f"{peak:<12}  (no matched data)")
            continue
        n       = len(sub)
        emp_m   = sub["empower_area"].mean()
        mc_m    = sub["area_cons"].mean()
        ratio   = mc_m / emp_m if emp_m > 0 else np.nan
        diff    = sub["area_diff_pct"].mean()
        u_pct   = (sub["U_area"] / sub["area_cons"] * 100).mean()
        cv      = sub["area_cv_pct"].mean()
        print(f"{peak:<12} {n:>4}  {emp_m:>14.1f}  {mc_m:>14.1f}  "
              f"{ratio:>7.4f}  {diff:>8.3f}%  {u_pct:>8.3f}%  {cv:>7.3f}%")

    print()

    # ── Summary table per peak — AIC model selection ──────────────────────────
    print("═══ mcescher vs Empower — area_AIC (AICc model selection) ═══\n")
    print(f"{'Peak':<12} {'n':>4}  {'Empower mean':>14}  "
          f"{'AIC mean':>14}  {'Ratio':>7}  {'Diff%':>8}  {'winner':>7}")
    print("─" * 80)

    for peak in NAMED_PEAKS:
        sub = merged[merged["name"] == peak].dropna(subset=["area_AIC", "empower_area"])
        if sub.empty:
            print(f"{peak:<12}  (no matched data)")
            continue
        n      = len(sub)
        emp_m  = sub["empower_area"].mean()
        aic_m  = sub["area_AIC"].mean()
        ratio  = aic_m / emp_m if emp_m > 0 else np.nan
        diff   = sub["aic_diff_pct"].mean()
        winner = sub["aic_winner"].mode()[0] if "aic_winner" in sub.columns else "?"
        print(f"{peak:<12} {n:>4}  {emp_m:>14.1f}  {aic_m:>14.1f}  "
              f"{ratio:>7.4f}  {diff:>8.3f}%  {winner:>7}")

    print()

    # ── RT comparison ─────────────────────────────────────────────────────────
    print(f"{'Peak':<12} {'RT Empower (min)':>18}  {'RT mcescher (min)':>18}  {'ΔRT (s)':>9}")
    print("─" * 65)
    for peak in NAMED_PEAKS:
        sub = merged[merged["name"] == peak].dropna(subset=["RT", "empower_RT_min"])
        if sub.empty:
            continue
        rt_mc  = sub["RT"].mean()
        rt_emp = sub["empower_RT_min"].mean()
        drt    = (rt_mc - rt_emp) * 60
        print(f"{peak:<12} {rt_emp:>18.4f}  {rt_mc:>18.4f}  {drt:>9.3f}")

    print()

    # ── Save merged table ─────────────────────────────────────────────────────
    out = os.path.join(data_dir, "comparison_table.csv")
    merged.to_csv(out, index=False)
    print(f"Full comparison table saved → {out}")


if __name__ == "__main__":
    data_dir = sys.argv[1] if len(sys.argv) > 1 else \
               "data/UNIS_FID_9000002_908779"
    main(data_dir)
