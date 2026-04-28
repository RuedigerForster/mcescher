"""
db/replicate_writer.py

Write one or all replicates from a processed dataset into MongoDB.

Collections written:
    samples     – one document per physical sample (upserted on sample_name + sample_type)
    replicates  – one document per CDF / injection

Usage:
    from db.replicate_writer import write_dataset
    from pymongo import MongoClient

    db = MongoClient("mongodb://localhost:27017")["mcescher"]
    write_dataset(
        data_dir   = "data/UNIS_FID_9000002_908779/",
        peaks_csv  = "data/UNIS_FID_9000002_908779/mcescher_peaks.csv",
        db         = db,
        method_ref = None,   # pass ObjectId once method collection is populated
    )
"""

import os
import re
from datetime import datetime, timedelta, timezone

import netCDF4 as nc
import numpy as np
import pandas as pd
from bson import ObjectId
from pymongo import MongoClient, ReturnDocument

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SENTINEL = 1e30

EMPOWER_NAME_MAP = {
    "CH4":       "Methane",
    "C2H6":      "Ethane",
    "C3H8":      "Propane",
    "i-C4H10":   "iso-Butane",
    "i-C4H10z":  "iso-Butane",
    "n-C4H10":   "n-Butane",
    "n-C4H10z":  "n-Butane",
}


# ---------------------------------------------------------------------------
# CDF helpers
# ---------------------------------------------------------------------------

def _parse_timestamp(stamp: str) -> datetime:
    """Parse Empower timestamp '20211230120342+0100' to UTC datetime."""
    dt_local = datetime.strptime(stamp[:14], "%Y%m%d%H%M%S")
    sign     = 1 if stamp[14] == "+" else -1
    offset   = timedelta(hours=int(stamp[15:17]), minutes=int(stamp[17:19]))
    return (dt_local - sign * offset).replace(tzinfo=timezone.utc)


def _safe_float(arr, i):
    v = float(arr[i])
    return None if v > SENTINEL else v


def read_cdf(cdf_path: str) -> dict:
    """Return metadata dict and Empower peak list from a CDF file."""
    ds = nc.Dataset(cdf_path)

    meta = {
        "acquired_at":        _parse_timestamp(ds.getncattr("injection_date_time_stamp")),
        "sample_name":        str(ds.getncattr("sample_name")).strip(),
        "sample_type":        str(ds.getncattr("sample_type")).strip(),
        "sample_id_comments": str(ds.getncattr("sample_id_comments")).strip(),
        "detector_unit":      str(ds.getncattr("detector_unit")).strip(),
        "detector_name":      str(ds.getncattr("detector_name")).strip(),
    }

    rt     = np.array(ds.variables["peak_retention_time"][:])
    area   = np.array(ds.variables["peak_area"][:])
    height = np.array(ds.variables["peak_height"][:])
    amount = np.array(ds.variables["peak_amount"][:])

    raw_names = ds.variables["peak_name"][:]
    names = []
    for row in raw_names:
        s = "".join(b.decode("latin-1") for b in row if b != b"\x00")
        names.append(re.sub(r"[^\x20-\x7E].*", "", s).strip())

    ds.close()

    empower_peaks = [
        {
            "empower_name":  names[i],
            "mapped_name":   EMPOWER_NAME_MAP.get(names[i]),
            "RT_s":          float(rt[i]),
            "RT_min":        float(rt[i]) / 60.0,
            "area":          _safe_float(area,   i),
            "height":        _safe_float(height, i),
            "amount":        _safe_float(amount, i),
        }
        for i in range(len(rt))
        if 0 < rt[i] < SENTINEL
    ]

    meta["empower_peaks"] = empower_peaks
    return meta


# ---------------------------------------------------------------------------
# Document builders
# ---------------------------------------------------------------------------

def _f(row, col):
    """Safe float from pandas row; returns None for NaN / missing."""
    v = row.get(col)
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)


def _s(row, col):
    """Safe string from pandas row."""
    v = row.get(col)
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else str(v)


def _empower_only_peak(ep: dict) -> dict:
    """Build a peak_table entry from Empower data alone (no mcescher)."""
    return {
        "RT":     ep["RT_min"],
        "name":   ep["mapped_name"],
        "height": ep["height"],
        "areas": {
            "PD": None, "TS": None, "gauss": None,
            "EGH": None, "AIC": None, "aic_winner": None, "cons": None,
        },
        "uncertainty": {
            "u_cons": None, "U_cons": None, "u_noise": None, "u_sass": None,
            "u_area": None, "U_area": None, "area_cv_pct": None,
            "Z_spread": None, "penalty": None,
        },
        "empower": {
            "RT_min": ep["RT_min"],
            "area":   ep["area"],
            "height": ep["height"],
        },
    }


def build_peak_table(mc_df: pd.DataFrame | None, empower_peaks: list) -> list:
    """
    Merge mcescher peaks with Empower peaks for one chromatogram.
    If mc_df is None (unprocessed file), returns Empower-only entries for
    named peaks.
    """
    emp_by_name = {
        ep["mapped_name"]: ep
        for ep in empower_peaks
        if ep["mapped_name"]
    }

    if mc_df is None or mc_df.empty:
        return [_empower_only_peak(ep) for ep in empower_peaks if ep["mapped_name"]]

    records = []
    for _, row in mc_df.iterrows():
        name = _s(row, "name")
        rec = {
            "RT":     _f(row, "RT"),
            "name":   name,
            "height": _f(row, "height"),
            "areas": {
                "PD":         _f(row, "area_PD"),
                "TS":         _f(row, "area_TS"),
                "gauss":      _f(row, "area_gauss"),
                "EGH":        _f(row, "area_EGH"),
                "AIC":        _f(row, "area_AIC"),
                "aic_winner": _s(row, "aic_winner"),
                "cons":       _f(row, "area_cons"),
            },
            "uncertainty": {
                "u_cons":      _f(row, "u_cons"),
                "U_cons":      _f(row, "U_cons"),
                "u_noise":     _f(row, "u_noise"),
                "u_sass":      _f(row, "u_sass"),
                "u_area":      _f(row, "u_area"),
                "U_area":      _f(row, "U_area"),
                "area_cv_pct": _f(row, "area_cv_pct"),
                "Z_spread":    _f(row, "Z_spread"),
                "penalty":     _f(row, "penalty"),
            },
        }
        if name and name in emp_by_name:
            ep = emp_by_name[name]
            rec["empower"] = {
                "RT_min": ep["RT_min"],
                "area":   ep["area"],
                "height": ep["height"],
            }
        records.append(rec)

    return records


def build_results(peak_table: list) -> list:
    """Named peaks only, initialised as accepted / not yet reported."""
    return [
        {
            "peak_name":      pk["name"],
            "aligned_RT":     pk["RT"],
            "reported_area":  pk["areas"]["cons"],
            "corrected_area": None,
            "areas":          pk["areas"],
            "uncertainty":    pk["uncertainty"],
            "empower":        pk.get("empower"),
            "status":         "accepted",
            "reported":       False,
        }
        for pk in peak_table
        if pk.get("name")
    ]


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

def write_replicate(cdf_path: str, db,
                    peaks_csv: str | None = None,
                    method_ref=None) -> ObjectId:
    """
    Insert one replicate document and upsert its parent sample.
    peaks_csv is optional; if absent, Empower data only is stored.
    Returns the inserted replicate _id.
    """
    source_file = os.path.basename(cdf_path)
    meta        = read_cdf(cdf_path)

    mc_df = None
    if peaks_csv and os.path.exists(peaks_csv):
        mc_all = pd.read_csv(peaks_csv)
        rows   = mc_all[mc_all["chrom"] == source_file]
        mc_df  = rows.copy() if not rows.empty else None

    peak_table = build_peak_table(mc_df, meta["empower_peaks"])
    results    = build_results(peak_table)

    # Upsert sample ─────────────────────────────────────────────────────────
    sample_doc = db.samples.find_one_and_update(
        filter={"sample_name": meta["sample_name"],
                "sample_type": meta["sample_type"]},
        update={"$setOnInsert": {
            "sample_name": meta["sample_name"],
            "sample_type": meta["sample_type"],
            "replicates":  [],
        }},
        upsert=True,
        return_document=ReturnDocument.AFTER,
    )

    # Insert replicate ───────────────────────────────────────────────────────
    rep_doc = {
        "sample_ref":         sample_doc["_id"],
        "method_ref":         method_ref,
        "source_file":        source_file,
        "acquired_at":        meta["acquired_at"],
        "sample_id_comments": meta["sample_id_comments"],
        "detector_name":      meta["detector_name"],
        "detector_unit":      meta["detector_unit"],
        "peak_table":         peak_table,
        "results":            results,
    }
    rep_id = db.replicates.insert_one(rep_doc).inserted_id

    # Push replicate key into sample ─────────────────────────────────────────
    db.samples.update_one(
        {"_id": sample_doc["_id"]},
        {"$push": {"replicates": rep_id}},
    )

    return rep_id


def write_dataset(data_dir: str, db,
                  peaks_csv: str | None = None,
                  method_ref=None) -> list[ObjectId]:
    """
    Write all CDF files found recursively under data_dir.
    peaks_csv is optional; matched by basename against chrom column.
    Returns list of inserted replicate _ids.
    """
    cdf_files = sorted(
        os.path.join(root, f)
        for root, _, files in os.walk(data_dir)
        for f in files
        if f.lower().endswith(".cdf")
    )

    inserted = []
    for cdf_path in cdf_files:
        rep_id = write_replicate(cdf_path, db, peaks_csv, method_ref)
        print(f"  {os.path.basename(cdf_path)}  →  {rep_id}")
        inserted.append(rep_id)

    print(f"\nInserted {len(inserted)} replicates.")
    return inserted


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Write mcescher replicates to MongoDB.")
    parser.add_argument("data_dir",          help="Root directory to scan for CDF files (recursive)")
    parser.add_argument("--peaks-csv",       default=None,
                        help="mcescher_peaks.csv (optional; Empower-only if omitted)")
    parser.add_argument("--uri",             default="mongodb://localhost:27017",
                        help="MongoDB connection URI")
    parser.add_argument("--db",              default="mcescher", help="Database name")
    args = parser.parse_args()

    client = MongoClient(args.uri)
    write_dataset(args.data_dir, client[args.db], peaks_csv=args.peaks_csv)
    client.close()
