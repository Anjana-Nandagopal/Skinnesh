# Skinnesh

**A reproducible tool to screen and rank plant compounds for skin-relevant activity and topical suitability.**

Skinnesh integrates structure standardisation, physicochemical screening, multi-target similarity scoring, ranking, leakage-controlled validation and HTML reporting in one command-line workflow.

## What makes this different

The proposal identified a gap between separate resources for plant compounds, bioactivity and property calculation. Skinnesh integrates those functions into a batch-capable workflow and keeps human and mushroom tyrosinase separate.

## Targets

- Human tyrosinase — ChEMBL CHEMBL1973
- Mushroom tyrosinase — ChEMBL CHEMBL3318 (*Agaricus bisporus*)
- Collagenase — ChEMBL CHEMBL333 (72 kDa type IV collagenase / MMP-2)
- Elastase — ChEMBL CHEMBL248 (human neutrophil elastase)

## Install

```bash
conda env create -f environment.yml
conda activate skinnesh
python -m skinnesh.cli --help
```

## 1. Build ChEMBL reference and validation sets

Internet access is required for this step.

```bash
python scripts/download_chembl.py --outdir data/chembl
```

The script retrieves public ChEMBL activity records, retains IC50/Ki measurements with usable numerical potency, labels pChEMBL >= 5 as active and <= 4 as inactive, and performs a scaffold-aware reference/validation split. The intermediate and summary files are written to `data/chembl/`.

**Important:** inspect `data/chembl/curation_summary.csv` before proceeding. If a target has insufficient data after curation, do not invent or silently substitute results.

## 2. Prepare the IMPPAT plant library

Download the chemical-information or 2D-structure file from the current IMPPAT release under its terms. Then identify its name and SMILES columns and run:

```bash
python scripts/prepare_imppat.py \
  --input <YOUR_IMPPAT_FILE> \
  --output data/plant_compounds.csv \
  --name-column <NAME_COLUMN> \
  --smiles-column <SMILES_COLUMN>
```

This script does not redistribute IMPPAT data; it converts the file you downloaded locally into the input schema used by Skinnesh.

## 3. Run the final pipeline

```bash
python scripts/run_final.py
```

This produces:

- `results/validation/validation_scores.csv`
- `results/validation/validation_metrics.csv`
- ROC and PR plots
- `results/final/ranked_compounds.csv`
- `results/final/report.html`

Open the final report:

```bash
open results/final/report.html
```

## Demo dataset

The small files in `data/` that were originally supplied with the project are for software testing only. **Do not use demo results as competition evidence.**

## Scoring

For each target, activity is the maximum Morgan-fingerprint/Tanimoto similarity to the target-specific reference-active set. The default combined score is:

`overall = 0.60 × mean target activity + 0.40 × topical suitability`

Weights are configurable in `skinnesh/config.yaml`.

The topical component is a transparent physicochemical heuristic using molecular weight, logP, TPSA, HBD and HBA. It is **not** a validated skin-permeation model.

## Validation

Validation is independent of the reference-active set and uses a scaffold-aware split. The pipeline reports ROC-AUC, PR-AUC and enrichment at the top 10% where both classes are available.

## Interpretation

A high Skinnesh score is a computational prioritisation signal. It does not prove enzyme inhibition, skin penetration, safety, efficacy or clinical/cosmetic benefit. Final biological interpretation must be reviewed by the participant.

## Final deliverables

See:

- `docs/REPORT_TEMPLATE.md`
- `docs/VIDEO_SCRIPT.md`
- `docs/SUBMISSION_CHECKLIST.md`
- `docs/DATA_SOURCES.md`
- `docs/AI_USAGE.md`

## AI disclosure

AI tools were used for brainstorming, code drafting, debugging assistance and documentation. Scientific decisions, data curation, validation design, biological interpretation and final conclusions were reviewed and determined by the participant.
