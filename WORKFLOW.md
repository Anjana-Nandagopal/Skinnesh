# Skinnesh workflow

1. **Public data acquisition**
   - ChEMBL target-specific enzyme activity records
   - IMPPAT phytochemical structures
2. **Activity curation**
   - IC50/Ki records with usable numerical potency
   - active: pChEMBL >= 5
   - inactive: pChEMBL <= 4
   - intermediate records excluded
3. **Leakage control**
   - scaffold-aware split into reference and held-out validation molecules
4. **Structure standardisation**
   - RDKit parsing, fragment-parent selection, uncharging, canonical SMILES
5. **Topical suitability**
   - MW, logP, TPSA, HBD, HBA heuristic score
6. **Target activity**
   - Morgan fingerprints + maximum Tanimoto similarity to target-specific reference actives
7. **Ranking**
   - configurable weighted combination of mean target activity and topical suitability
8. **Validation**
   - ROC-AUC, PR-AUC, enrichment factor at 10%
9. **Reporting**
   - ranked CSV, validation CSV/plots, HTML report
10. **Biological interpretation**
   - participant-reviewed interpretation of the final candidates and target patterns
