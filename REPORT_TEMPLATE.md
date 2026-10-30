# Skinnesh: Screening and ranking plant compounds for skin-relevant activity and topical suitability

## Abstract

Skinnesh is a reproducible Python workflow designed to prioritise plant-derived compounds for activity across human tyrosinase, mushroom tyrosinase, collagenase and elastase while incorporating a transparent topical-suitability heuristic. The workflow standardises input structures, calculates physicochemical descriptors, computes Morgan-fingerprint similarity to curated reference actives, combines activity and topical scores, and evaluates the similarity-ranking strategy using a held-out validation set.

**Final result summary:** [INSERT ACTUAL NUMBER OF PLANT COMPOUNDS, VALIDATION COUNTS, AUCs, ENRICHMENT, AND TOP CANDIDATES.]

## 1. Problem and research gap

Plant-derived compounds are widely investigated for pigmentation and skin-ageing related applications. Existing resources provide plant-compound structures, experimental bioactivity and physicochemical-property calculations, but these components are generally separate. The proposed gap is an integrated, batch-capable and reproducible workflow that combines multi-target activity prioritisation, topical-suitability scoring and independent validation while keeping human and mushroom tyrosinase distinct.

## 2. Methods

### 2.1 Input standardisation

SMILES are parsed with RDKit, fragment parents are selected, charges are neutralised where possible, and canonical SMILES are generated. Invalid structures are excluded from downstream scoring and reported.

### 2.2 Physicochemical descriptors and topical suitability

Skinnesh calculates molecular weight, logP, TPSA, hydrogen-bond donors, hydrogen-bond acceptors, rotatable bonds and aromatic rings. A transparent heuristic combines selected descriptor bands into a 0–1 topical-suitability score. This score is a screening heuristic and is not a validated skin-permeation model.

### 2.3 Target activity scoring

For each target, Morgan fingerprints are generated and the maximum Tanimoto similarity to the target-specific reference-active set is used as the activity-prioritisation score.

### 2.4 Combined ranking

The current default is:

`overall = 0.60 × mean target-activity score + 0.40 × topical-suitability score`

The weights are configurable.

### 2.5 Data curation and leakage control

ChEMBL records are filtered to IC50 and Ki measurements with usable numerical activity information. The default activity labels use pChEMBL ≥ 5 as active and pChEMBL ≤ 4 as inactive, excluding the intermediate region. Reference and validation molecules are separated by molecular scaffold to reduce analogue leakage.

### 2.6 Validation

Validation is performed independently for each target. ROC-AUC, PR-AUC and enrichment factor at the top 10% are reported where both classes are available.

## 3. Results

### 3.1 Dataset summary

[INSERT `data/chembl/curation_summary.csv` TABLE.]

### 3.2 Validation

[INSERT `results/validation/validation_metrics.csv` TABLE AND ROC/PR FIGURES.]

### 3.3 Plant-compound ranking

[INSERT TOP 10–20 TABLE FROM `results/final/ranked_compounds.csv`.]

### 3.4 Biological interpretation

Write the interpretation yourself after reviewing the final ranked molecules and the literature. Discuss target-specific patterns, multi-target candidates and any important limitations. Do not treat computational similarity as experimental proof.

## 4. Limitations

- Similarity to known actives is not a mechanistic proof of inhibition.
- The topical-suitability component is a transparent heuristic, not a skin-permeation or formulation model.
- ChEMBL assays differ in experimental design, substrate and conditions; target-level curation reduces but does not eliminate heterogeneity.
- Validation performance depends on the curated reference and held-out sets.
- Plant origin does not imply safety, efficacy or suitability for cosmetic use.

## 5. Reproducibility

The repository contains an environment file, command-line interface, configuration, data-preparation scripts, validation workflow and HTML reporting. A fresh environment should reproduce the computational outputs when the same public input files and versions are used.

## 6. AI usage declaration

AI tools were used for brainstorming, code drafting, debugging assistance and documentation. Scientific decisions, dataset curation, validation design, interpretation and final conclusions were reviewed and determined by the participant.

## 7. References and data sources

- ChEMBL — public bioactivity database; record target IDs and retrieval date in the final repository.
- IMPPAT — Indian Medicinal Plants, Phytochemistry And Therapeutics; record the exact release/file and applicable licence.
- RDKit — cheminformatics toolkit.
- Scikit-learn — validation metrics.
