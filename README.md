# GeneSelectR 2.0: benchmark and supplementary material

This repository accompanies the manuscript **"GeneSelectR 2.0: An R Workflow
for Predictive Gene Selection and Biological Interpretation"**. It contains the
source tables behind every manuscript table and figure, the scripts that
produced them, and the supplementary analyses referred to in the text. The
GeneSelectR 2.0 R package itself is provided in a separate repository.

## Manuscript to repository map

| Manuscript item | Source tables (`manuscript/source_data/`) | Produced by |
|---|---|---|
| Table 1: datasets | — | `benchmarks/validation_datasets.R`, `benchmarks/validation_prepare.R`, `redesign/R/*_data.R` |
| Figure 1: workflow | — | `manuscript/make_workflow_figure.py` |
| Table 2, Figure 2: predictive performance | `by_dataset_primary.csv`, `by_dataset_size.csv`, `by_validation_*.csv`, `random30_by_dataset_size.csv`, `prediction_table2_verified.csv` | `supplementary/scripts/run_primary.sh`, `run_random_reference_sensitivity.sh`; Stabl: `redesign/stabl_comparison/run_comparison.R` |
| Figure 3: selection stability | `stability_classical.csv`, `stability_by_dataset_primary.csv` | `supplementary/scripts/run_stability.sh`; Stabl: `redesign/stabl_comparison/summarise.R` |
| Section 4.3: biological enrichment | `biology_comparison_by_dataset.csv`, `biology_by_dataset.csv`, `ot_*.csv`, `stabl_biology_*.csv`, `biology_enrichment_counts.csv` | `supplementary/scripts/run_biology.sh`; Stabl: `redesign/stabl_comparison/run_biology.R` |
| Section 4.4: training-based size selection | reported in [`supplementary/results_summary.md`](supplementary/results_summary.md) | `supplementary/scripts/run_gene_set_size_sensitivity.sh` |
| Section 4.5: component analysis | reported in [`supplementary/synthesized_report.md`](supplementary/synthesized_report.md) | `supplementary/scripts/run_component_ablation.sh`, `run_permutation_adjustment.sh` |
| Figure 4: asthma case study | `asthma_eight.csv`, `asthma_biology.csv`, `asthma_k20_biology_comparison.csv`, `by_dataset_size.csv` | `redesign/complementarity_analysis/03_gene_group_characterization.R`, `supplementary/scripts/run_biology.sh` |

## Repository layout

```text
├── manuscript/                 Source tables and scripts for the manuscript figures and tables
│   ├── source_data/            Aggregated benchmark and biological results (CSV)
│   ├── figures/                Figures 1–4 (PDF, PNG)
│   ├── make_figures.py         Figures 2–4 from source_data/
│   └── make_workflow_figure.py Figure 1
├── supplementary/              Supplementary material
│   ├── methods/                Supplementary methods
│   ├── results_summary.md      Numerical results, including ablations
│   ├── synthesized_report.md   Integrated analysis report
│   ├── figures/                Supplementary figures S1–S6
│   ├── configurations/         Registry of reported, sensitivity and development configurations
│   ├── provenance/             Revision and reproducibility records
│   └── scripts/                Ordered entry points for the full re-analysis
├── redesign/                   Analysis code
│   ├── R/                      Shared fitting, evaluation and biological functions
│   ├── stabl_comparison/       Stabl benchmark and biological assessment
│   ├── complementarity_analysis/  Stability, DGE comparison and gene-level analyses
│   ├── run_older7_*.R          Ablation, sensitivity and biology runners used in the paper
│   └── run_*.R, tests/         Primary runners and method-development experiments
├── benchmarks/                 Dataset preparation and earlier benchmark implementations
└── analysis/                   Shared configuration, input checks and package bootstrap
```

Scripts under `redesign/` that are not called from `supplementary/scripts/`
record method development (alternative priors, score combinations, pruning,
greedy panels and similar experiments). They are kept for transparency and are
not needed to reproduce the manuscript. The
[configuration registry](supplementary/configurations/registry.md) classifies
each one.

## Reproducing the manuscript figures

The manuscript figures can be regenerated from the committed source tables
without the raw data:

```bash
cd manuscript
python3 make_figures.py          # Figures 2–4 (NumPy, pandas, Matplotlib)
python3 make_workflow_figure.py  # Figure 1 (NumPy, Matplotlib, cairosvg)
```

## Full re-analysis

Run from the repository root, in this order:

```bash
bash supplementary/scripts/check_inputs.sh      # lists missing inputs and packages
bash supplementary/scripts/run_tests.sh
bash supplementary/scripts/run_primary.sh       # seven-dataset benchmark
bash supplementary/scripts/run_biology.sh
bash supplementary/scripts/run_stability.sh
bash supplementary/scripts/run_component_ablation.sh
bash supplementary/scripts/run_permutation_adjustment.sh
bash supplementary/scripts/run_gene_set_size_sensitivity.sh
bash supplementary/scripts/run_random_reference_sensitivity.sh
bash supplementary/scripts/run_stabl_comparison.sh   # requires the Stabl Python environment
```

`bash supplementary/scripts/run_all.sh` runs the R-based scopes in sequence.
The complete analysis takes many hours. Two workers are used by default; set
`GENESELECTR_N_CORES` to change this. `bash analysis/bootstrap_package.sh`
extracts the pinned GeneSelectR 2.0 package revision (`621be0c1`) used by the
runners.

## Data

| Dataset | Source | Task |
|---|---|---|
| SOS-ALL (E-MTAB-16535) | ArrayExpress | Atopic dermatitis / healthy (development) |
| GSE101794 | GEO | Crohn's disease / non-IBD |
| GSE107994 | GEO | Active / latent tuberculosis |
| GSE13355 | GEO | Lesional / non-lesional psoriasis |
| GSE65682 | GEO | Death by day 28 / survival in sepsis |
| GSE69683 | GEO | Severe / moderate asthma |
| IMvigor210 | `easierData` (Bioconductor) | Immunotherapy response / non-response |

Expression data are not redistributed. The preparation scripts in
`benchmarks/` and `redesign/R/` download or read them from the public sources
into the ignored `data/` directory. Frozen annotation inputs for the Stabl
biological assessment are included under
`redesign/stabl_comparison/biology/reference_inputs/`.

## License

See [LICENSE](LICENSE).
