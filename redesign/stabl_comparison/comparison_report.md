# Empirical Stabl comparison

## Design and verification

Stabl was fitted independently in each of the 15 saved outer-training divisions of the six validation datasets and the SOS-ALL development dataset. The same 2,000-gene candidate pools, training-based standardization, panel sizes (10, 20, 50, 100, 200, and 500 genes), outer-test samples, ensemble evaluator, and evaluation seeds were used as in the existing benchmark. Existing GeneSelectR 2.0 and classical-method results were read from the saved benchmark files. No existing feature-selection method was refitted.

The analysis produced 105 Stabl rankings and 630 outer-test AUC values. Every ranking contained 2,000 distinct candidate genes; every panel contained the requested number of genes; all AUCs were finite. The evaluator reproduced two saved GeneSelectR 2.0 AUCs for GSE101794 division 1/1: 0.9745098 at 10 genes and 0.9843137 at 500 genes. The recalculated GeneSelectR 2.0 Nogueira stability values matched the manuscript figure to three decimal places in all seven datasets.

## Primary comparison

AUC is averaged over 15 outer-test divisions and gene-set sizes of 10, 20, and 50. Nogueira stability is averaged over the same three sizes. Each validation dataset contributes equally to the validation mean. ΔAUC is GeneSelectR 2.0 minus Stabl.

| Dataset | GeneSelectR AUC | Stabl AUC | ΔAUC | GeneSelectR stability | Stabl stability |
|---|---:|---:|---:|---:|---:|
| GSE101794, Crohn disease | 0.9792 | 0.9805 | −0.0014 | 0.605 | 0.607 |
| GSE107994, tuberculosis | 0.9679 | 0.9659 | +0.0020 | 0.616 | 0.628 |
| GSE13355, psoriasis | 0.9968 | 0.9982 | −0.0013 | 0.449 | 0.656 |
| GSE65682, sepsis | 0.6645 | 0.6704 | −0.0059 | 0.332 | 0.479 |
| GSE69683, asthma | 0.7872 | 0.7719 | +0.0153 | 0.526 | 0.468 |
| IMvigor210, treatment response | 0.6700 | 0.6176 | +0.0523 | 0.320 | 0.307 |
| **Six-validation-dataset mean** | **0.8443** | **0.8341** | **+0.0102** | **0.475** | **0.524** |
| SOS-ALL, development | 0.6344 | 0.6265 | +0.0079 | 0.305 | 0.340 |

Stabl had higher mean AUC in three of six validation datasets and higher stability in four. The six-dataset mean AUC difference was 0.0102 in favor of GeneSelectR 2.0. The mean Nogueira stability difference was 0.050 in favor of Stabl. These are descriptive comparisons across correlated repeated folds; no significance test was applied.

## All benchmark gene-set sizes

| Genes | GeneSelectR mean AUC | Stabl mean AUC | ΔAUC |
|---:|---:|---:|---:|
| 10 | 0.8366 | 0.8267 | +0.0099 |
| 20 | 0.8489 | 0.8376 | +0.0113 |
| 50 | 0.8472 | 0.8379 | +0.0093 |
| 100 | 0.8400 | 0.8404 | −0.0003 |
| 200 | 0.8411 | 0.8454 | −0.0044 |
| 500 | 0.8467 | 0.8478 | −0.0011 |

The fixed-size Stabl ranking uses maximum selection frequency across its regularization grid. Stabl's native data-derived threshold selected a median of 4 genes per outer division (range 0–111); that threshold was recorded separately. At 500 genes, 43 of 105 panels included genes with zero selection frequency, resolved by the saved training-variance order. All 10-, 20-, 50-, 100-, and 200-gene panels comprised genes with positive selection frequency. Results at 500 genes therefore have a larger contribution from the fixed-size extension.

## Interpretation and limits

The primary comparison evaluates Stabl's ranking at the same panel sizes as GeneSelectR 2.0. It does not evaluate Stabl's variable-size native output as a separate predictive strategy. The original and extension runs used the same evaluator-package versions; the original manifest records R 4.5.2 and the extension used R 4.6.1. Two final-pass `glmnet` warnings concerned late regularization-path convergence in SOS-ALL. All resulting probabilities and AUCs were finite.

MVFS-SHAP was not evaluated. No public authors' code or package was identified from the [publication record](https://pubmed.ncbi.nlm.nih.gov/41289809/) or related code search. The Stabl code and configuration are documented in [the authors' repository](https://github.com/gregbellan/Stabl).
