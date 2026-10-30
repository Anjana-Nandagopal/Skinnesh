from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

TARGETS = ["human_tyrosinase", "mushroom_tyrosinase", "collagenase", "elastase"]


def standardize_smiles(smiles: str) -> Tuple[str | None, str]:
    if pd.isna(smiles):
        return None, "missing"
    smiles = str(smiles).strip()
    if not smiles:
        return None, "empty"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "invalid_smiles"
    try:
        mol = rdMolStandardize.FragmentParent(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, canonical=True), "ok"
    except Exception:
        return None, "standardization_failed"


def descriptor_row(mol):
    return {
        "MW": Descriptors.MolWt(mol),
        "logP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "RotB": Lipinski.NumRotatableBonds(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    }


def _band_score(value, ideal_min=None, ideal_max=None, hard_min=None, hard_max=None):
    if value is None or not np.isfinite(value):
        return 0.0
    if ideal_min is not None and value < ideal_min:
        if hard_min is None:
            return 0.0
        return max(0.0, (value - hard_min) / (ideal_min - hard_min))
    if ideal_max is not None and value > ideal_max:
        if hard_max is None:
            return 0.0
        return max(0.0, (hard_max - value) / (hard_max - ideal_max))
    return 1.0


def topical_score(desc: Dict, cfg: Dict) -> float:
    t = cfg["topical"]
    parts = [
        _band_score(desc["MW"], ideal_max=t["mw"]["ideal_max"], hard_max=t["mw"]["hard_max"]),
        _band_score(desc["logP"], ideal_min=t["logp"]["ideal_min"], ideal_max=t["logp"]["ideal_max"], hard_min=t["logp"].get("hard_min", -1), hard_max=t["logp"].get("hard_max", 7)),
        _band_score(desc["TPSA"], ideal_min=t["tpsa"]["ideal_min"], ideal_max=t["tpsa"]["ideal_max"], hard_min=t["tpsa"].get("hard_min", 5), hard_max=t["tpsa"].get("hard_max", 150)),
        _band_score(desc["HBD"], ideal_max=t["hbd"]["ideal_max"], hard_max=t["hbd"].get("hard_max", 10)),
        _band_score(desc["HBA"], ideal_max=t["hba"]["ideal_max"], hard_max=t["hba"].get("hard_max", 20)),
    ]
    return float(np.mean(parts))


def morgan_fp(mol, radius=2, n_bits=2048):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return generator.GetFingerprint(mol)


def max_similarity(query_fp, ref_fps) -> float:
    if not ref_fps:
        return 0.0
    sims = DataStructs.BulkTanimotoSimilarity(query_fp, ref_fps)
    return float(max(sims)) if sims else 0.0


def load_reference_data(path: str | Path):
    df = pd.read_csv(path)
    required = {"target", "smiles"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Reference file missing columns: {sorted(missing)}")
    refs = {}
    for target, group in df.groupby("target"):
        fps = []
        for s in group["smiles"].dropna():
            canonical, _ = standardize_smiles(s)
            if canonical:
                mol = Chem.MolFromSmiles(canonical)
                if mol:
                    fps.append(morgan_fp(mol))
        refs[target] = fps
    return refs


def rank_compounds(input_path, reference_path, output_dir, cfg):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    compounds = pd.read_csv(input_path)
    if "smiles" not in compounds.columns:
        raise ValueError("Input CSV must contain a 'smiles' column.")
    if "name" not in compounds.columns:
        compounds.insert(0, "name", [f"compound_{i+1}" for i in range(len(compounds))])
    refs = load_reference_data(reference_path)
    rows = []
    for _, row in compounds.iterrows():
        canonical, status = standardize_smiles(row["smiles"])
        out = row.to_dict()
        out["standardization_status"] = status
        if canonical is None:
            out.update({"canonical_smiles": "", "MW": np.nan, "logP": np.nan, "TPSA": np.nan,
                        "HBD": np.nan, "HBA": np.nan, "RotB": np.nan, "AromaticRings": np.nan,
                        "topical_suitability": 0.0})
            for target in TARGETS:
                out[f"{target}_activity"] = 0.0
            out["activity_mean"] = 0.0
            out["overall_score"] = 0.0
            rows.append(out)
            continue
        mol = Chem.MolFromSmiles(canonical)
        desc = descriptor_row(mol)
        out["canonical_smiles"] = canonical
        out.update(desc)
        out["topical_suitability"] = topical_score(desc, cfg)
        qfp = morgan_fp(mol)
        activity_values = []
        for target in TARGETS:
            score = max_similarity(qfp, refs.get(target, []))
            out[f"{target}_activity"] = score
            activity_values.append(score)
        out["activity_mean"] = float(np.mean(activity_values))
        out["overall_score"] = (
            cfg["weights"]["activity"] * out["activity_mean"]
            + cfg["weights"]["topical"] * out["topical_suitability"]
        )
        rows.append(out)
    result = pd.DataFrame(rows).sort_values(["overall_score", "activity_mean"], ascending=False)
    result["rank"] = np.arange(1, len(result) + 1)
    result.to_csv(output_dir / "ranked_compounds.csv", index=False)
    return result
