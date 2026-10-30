from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve
from .core import load_reference_data, morgan_fp, max_similarity, standardize_smiles


def enrichment_at_k(y_true, scores, k_frac=0.10):
    n = len(y_true)
    if n == 0 or np.sum(y_true) == 0:
        return np.nan
    k = max(1, int(np.ceil(n * k_frac)))
    order = np.argsort(-np.asarray(scores))[:k]
    top_rate = np.mean(np.asarray(y_true)[order])
    base_rate = np.mean(y_true)
    return float(top_rate / base_rate) if base_rate else np.nan


def score_validation_set(validation_path, reference_path, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    val = pd.read_csv(validation_path)
    required = {"smiles", "target", "label"}
    missing = required - set(val.columns)
    if missing:
        raise ValueError(f"Validation file missing columns: {sorted(missing)}")
    refs = load_reference_data(reference_path)
    scored = []
    for _, row in val.iterrows():
        canonical, status = standardize_smiles(row["smiles"])
        if not canonical:
            continue
        mol = Chem.MolFromSmiles(canonical)
        scored.append({
            "name": row.get("name", ""),
            "target": row["target"],
            "label": int(row["label"]),
            "similarity_score": max_similarity(morgan_fp(mol), refs.get(row["target"], [])),
            "standardization_status": status,
        })
    scored = pd.DataFrame(scored)
    metrics = []
    for target, group in scored.groupby("target"):
        if group["label"].nunique() < 2:
            continue
        y = group["label"].to_numpy()
        s = group["similarity_score"].to_numpy()
        metrics.append({
            "target": target,
            "n": len(group),
            "actives": int(y.sum()),
            "decoys": int((y == 0).sum()),
            "roc_auc": roc_auc_score(y, s),
            "pr_auc": average_precision_score(y, s),
            "enrichment_factor_at_10pct": enrichment_at_k(y, s, 0.10),
        })
    metrics_df = pd.DataFrame(metrics)
    scored.to_csv(output_dir / "validation_scores.csv", index=False)
    metrics_df.to_csv(output_dir / "validation_metrics.csv", index=False)
    import matplotlib.pyplot as plt
    for target, group in scored.groupby("target"):
        if group["label"].nunique() < 2:
            continue
        y = group["label"].to_numpy(); s = group["similarity_score"].to_numpy()
        fpr, tpr, _ = roc_curve(y, s)
        auc = roc_auc_score(y, s)
        plt.figure(figsize=(6, 5)); plt.plot(fpr, tpr, label=f"ROC-AUC={auc:.3f}")
        plt.plot([0, 1], [0, 1], linestyle="--")
        plt.xlabel("False positive rate"); plt.ylabel("True positive rate")
        plt.title(f"Skinnesh validation — {target}"); plt.legend(); plt.tight_layout()
        plt.savefig(output_dir / f"roc_{target}.png", dpi=180); plt.close()
        precision, recall, _ = precision_recall_curve(y, s)
        plt.figure(figsize=(6, 5)); plt.plot(recall, precision)
        plt.xlabel("Recall"); plt.ylabel("Precision")
        plt.title(f"Skinnesh precision-recall — {target}"); plt.tight_layout()
        plt.savefig(output_dir / f"pr_{target}.png", dpi=180); plt.close()
    return metrics_df
