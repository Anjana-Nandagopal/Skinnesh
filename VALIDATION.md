# Validation design

Skinnesh separates the molecules used to construct the reference-active fingerprints from the molecules used to evaluate ranking performance.

The ChEMBL preparation script:

1. retrieves public target-specific activity records;
2. retains IC50 and Ki measurements with usable potency values;
3. uses pChEMBL where available and otherwise converts standard values to nM;
4. labels pChEMBL >= 5 as active and <= 4 as inactive;
5. removes the intermediate 4–5 region;
6. deduplicates by target and ChEMBL molecule;
7. computes Bemis–Murcko-style molecular scaffolds with RDKit;
8. assigns whole scaffolds to an approximately 80/20 reference/validation split.

The held-out set is scored only against reference actives. Metrics are computed separately by target:

- ROC-AUC
- PR-AUC
- enrichment factor at the top 10%

If a target does not have sufficient held-out examples in both classes, Skinnesh reports that rather than fabricating or silently substituting a metric.

## Important caveat

The active/inactive thresholds are operational screening definitions for this workflow, not universal biological definitions. The final report should describe them explicitly and discuss assay heterogeneity.
