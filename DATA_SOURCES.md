# Data sources

## ChEMBL

Used for target-specific experimental bioactivity records. The current target mapping in the pipeline is:

| Skinnesh channel | ChEMBL target | Description |
|---|---|---|
| human_tyrosinase | CHEMBL1973 | Human tyrosinase |
| mushroom_tyrosinase | CHEMBL3318 | Tyrosinase from *Agaricus bisporus* |
| collagenase | CHEMBL333 | 72 kDa type IV collagenase / MMP-2 |
| elastase | CHEMBL248 | Human neutrophil/neutrophil elastase |

Record the retrieval date in the final report because ChEMBL is updated over time.

## IMPPAT

The plant screening library is sourced from IMPPAT. Use the exact release and downloaded file in the final report. The 2026 IMPPAT 3.0 release provides batch downloads including chemical information and 2D structures. Check the current licence/terms before redistributing the raw dataset.

## Software

- RDKit for chemical standardisation, descriptors and fingerprints.
- pandas/numpy for data processing.
- scikit-learn for validation metrics.
- matplotlib for plots.
