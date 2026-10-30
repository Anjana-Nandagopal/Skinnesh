#!/usr/bin/env python
"""Robust ChEMBL downloader/curator for Skinnesh.

Downloads only IC50 and Ki activity records in small pages with retries,
curates usable structures/measurements, and creates scaffold-separated
reference-active and validation sets.
"""
from __future__ import annotations
import argparse, hashlib, time
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

TARGETS = {
    "human_tyrosinase": "CHEMBL1973",
    "mushroom_tyrosinase": "CHEMBL3318",
    "collagenase": "CHEMBL333",
    "elastase": "CHEMBL248",
}
API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"


def make_session():
    s = requests.Session()
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        status=5,
        backoff_factor=2,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        respect_retry_after_header=True,
    )
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.headers.update({"User-Agent": "Skinnesh-GAP2TOOL/1.0"})
    return s


SESSION = make_session()


def get_json(params, attempts=4):
    last = None
    for attempt in range(1, attempts + 1):
        try:
            r = SESSION.get(API, params=params, timeout=(20, 180))
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            last = e
            if attempt < attempts:
                wait = min(30, 2 ** attempt)
                print(f"  Request failed (attempt {attempt}/{attempts}): {e}")
                print(f"  Retrying in {wait}s...")
                time.sleep(wait)
    raise RuntimeError(f"ChEMBL request failed after {attempts} attempts: {last}")


def download_target(target_name, target_id, out_csv, page_size=100, max_rows=5000):
    rows = []
    # Query IC50 and Ki separately so the server returns a much smaller result set.
    for standard_type in ["IC50", "Ki"]:
        offset = 0
        print(f"  Downloading {target_name} {standard_type} records...")
        while len(rows) < max_rows:
            remaining = max_rows - len(rows)
            limit = min(page_size, remaining)
            params = {
                "target_chembl_id": target_id,
                "standard_type": standard_type,
                "limit": limit,
                "offset": offset,
            }
            data = get_json(params)
            batch = data.get("activities", [])
            if not batch:
                break
            rows.extend(batch)
            offset += len(batch)
            print(f"    fetched {len(rows)} records so far")
            page_meta = data.get("page_meta", {})
            if len(batch) < limit or not page_meta.get("next"):
                break
            time.sleep(0.2)
        if len(rows) >= max_rows:
            break

    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError(f"No activity records returned for {target_name} ({target_id})")
    df["skinnesh_target"] = target_name
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df


def standardize(smiles):
    if pd.isna(smiles) or str(smiles).strip().lower() in {"", "nan", "none"}:
        return None
    try:
        m = Chem.MolFromSmiles(str(smiles))
        if m is None:
            return None
        return Chem.MolToSmiles(m, canonical=True)
    except Exception:
        return None


def unit_to_nM(value, unit):
    if pd.isna(value) or pd.isna(unit):
        return np.nan
    try:
        v = float(value)
    except Exception:
        return np.nan
    u = str(unit).strip().lower().replace("µ", "u").replace("μ", "u")
    factors = {
        "m": 1e9, "mm": 1e6, "um": 1e3, "nm": 1.0,
        "pm": 1e-3, "mol/l": 1e9, "mol/liter": 1e9,
    }
    return v * factors.get(u, np.nan)


def curate(raw):
    df = raw.copy()
    df = df[df["standard_type"].isin(["IC50", "Ki"])].copy()
    # Keep exact numeric measurements only; qualifiers such as > or < are excluded.
    if "standard_relation" in df.columns:
        df = df[df["standard_relation"].fillna("=").eq("=")].copy()

    df["smiles"] = df.get("canonical_smiles", pd.Series(index=df.index, dtype=object))
    df["smiles"] = df["smiles"].map(standardize)
    df = df[df["smiles"].notna()].copy()

    df["pchembl_value"] = pd.to_numeric(df.get("pchembl_value"), errors="coerce")
    df["value_nM"] = [
        unit_to_nM(v, u)
        for v, u in zip(df.get("standard_value"), df.get("standard_units"))
    ]

    calculated = pd.Series(np.nan, index=df.index, dtype=float)
    valid_nM = df["value_nM"].notna() & (df["value_nM"] > 0)
    calculated.loc[valid_nM] = -np.log10(df.loc[valid_nM, "value_nM"] * 1e-9)
    df["pchembl_calc"] = df["pchembl_value"].where(df["pchembl_value"].notna(), calculated)
    df = df[np.isfinite(df["pchembl_calc"])].copy()

    # Conservative binary labels: <=10 uM active; >=100 uM inactive.
    df["label"] = df["pchembl_calc"].apply(
        lambda x: 1 if x >= 5.0 else (0 if x <= 4.0 else -1)
    )
    df = df[df["label"] >= 0].copy()
    df["molecule_chembl_id"] = df.get("molecule_chembl_id", "")
    df = df.drop_duplicates(
        subset=["skinnesh_target", "molecule_chembl_id", "standard_type"],
        keep="first",
    )
    keep = [
        "skinnesh_target", "molecule_chembl_id", "smiles", "standard_type",
        "standard_relation", "standard_value", "standard_units", "pchembl_calc", "label",
    ]
    return df[[c for c in keep if c in df.columns]]


def scaffold_key(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return ""
    scaf = MurckoScaffold.GetScaffoldForMol(m)
    return Chem.MolToSmiles(scaf, canonical=True) if scaf.GetNumAtoms() > 0 else Chem.MolToSmiles(m, canonical=True)


def scaffold_split(df, seed="skinnesh2026"):
    df = df.copy()
    df["scaffold"] = df["smiles"].map(scaffold_key)
    keys = sorted(df["scaffold"].unique(), key=lambda x: hashlib.sha256((seed + x).encode()).hexdigest())
    if len(keys) <= 1:
        return df.copy(), df.iloc[0:0].copy()
    n_train = max(1, int(round(len(keys) * 0.8)))
    n_train = min(n_train, len(keys) - 1)
    train_scaf = set(keys[:n_train])
    train = df[df.scaffold.isin(train_scaf)].copy()
    test = df[~df.scaffold.isin(train_scaf)].copy()
    return train, test


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/chembl")
    ap.add_argument("--max-rows", type=int, default=5000)
    ap.add_argument("--min-per-class", type=int, default=10)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    all_curated = []
    summaries = []

    for name, tid in TARGETS.items():
        print(f"\n=== {name} ({tid}) ===")
        raw_path = out / f"{name}_raw.csv"
        try:
            raw = download_target(name, tid, raw_path, max_rows=args.max_rows)
            cur = curate(raw)
            train, test = scaffold_split(cur)
            ref = train[train.label == 1].copy()
            val = test.copy()

            status = "OK"
            if len(ref) < args.min_per_class or val.label.nunique() < 2:
                status = "INSUFFICIENT_AFTER_SPLIT"
            else:
                ref[["skinnesh_target", "smiles"]].rename(
                    columns={"skinnesh_target": "target"}
                ).to_csv(out / f"{name}_reference_actives.csv", index=False)
                val[["skinnesh_target", "smiles", "label"]].rename(
                    columns={"skinnesh_target": "target"}
                ).assign(name=val["molecule_chembl_id"].astype(str)).to_csv(
                    out / f"{name}_validation.csv", index=False
                )
                all_curated.append(cur)

            summaries.append({
                "target": name,
                "raw": len(raw),
                "curated": len(cur),
                "reference_actives": len(ref),
                "validation_n": len(val),
                "validation_actives": int((val.label == 1).sum()),
                "validation_inactives": int((val.label == 0).sum()),
                "status": status,
            })
        except Exception as e:
            print(f"ERROR for {name}: {e}")
            summaries.append({"target": name, "raw": 0, "curated": 0, "reference_actives": 0,
                              "validation_n": 0, "validation_actives": 0,
                              "validation_inactives": 0, "status": f"DOWNLOAD_ERROR: {e}"})

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(out / "curation_summary.csv", index=False)
    print("\n=== CURATION SUMMARY ===")
    print(summary_df.to_string(index=False))

    good = summary_df["status"].eq("OK").all() and len(summary_df) == len(TARGETS)
    if good:
        refs = pd.concat([pd.read_csv(out / f"{n}_reference_actives.csv") for n in TARGETS], ignore_index=True)
        vals = pd.concat([pd.read_csv(out / f"{n}_validation.csv") for n in TARGETS], ignore_index=True)
        refs.to_csv("data/reference_actives.csv", index=False)
        vals.to_csv("data/validation.csv", index=False)
        print("\nCreated data/reference_actives.csv and data/validation.csv")
    else:
        print("\nSome targets did not produce sufficient data; no combined final files were created.")


if __name__ == "__main__":
    main()
