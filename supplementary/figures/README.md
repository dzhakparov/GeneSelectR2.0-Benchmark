# Supplementary figures

## Verification state

Six rendered figures are committed in this directory. They are the archived
renderings produced by the maintained
`redesign/complementarity_analysis/06_make_figures.R` workflow. The figure
source tables contain the seven primary datasets and the four follow-up
datasets retained in the earlier figure set. Source CSV tables remain in the
ignored result directory because generated tabular outputs are excluded from
the repository commit.

## Figure inventory

| Figure | Files | Generation script | Caption |
|---|---|---|---|
| Figure 2: prediction and stability | `Figure_2_prediction_and_stability.pdf` / `.png` | `redesign/complementarity_analysis/06_make_figures.R` | Mean test AUC and cross-split gene-set stability by method, dataset, and panel size |
| Supplementary stability curves and heatmap | `Supplementary_stability_by_gene_set_size.pdf` / `.png`; `Supplementary_stability_heatmap_k20.pdf` / `.png` | `redesign/complementarity_analysis/06_make_figures.R` | Stability measures across panel sizes and the archived dataset set |
| Figure 3: DGE/GeneSelectR complementarity | `Figure_3_DGE_GeneSelectR_complementarity.pdf` / `.png` | `redesign/complementarity_analysis/06_make_figures.R` | Rank percentiles, top-k overlap, and recurrent genes from matched candidate pools |
| Supplementary DGE score components | `Supplementary_DGE_vs_GeneSelectR_components.pdf` / `.png` | `redesign/complementarity_analysis/06_make_figures.R` | DGE rank compared with GeneSelectR score components |
| Figure 4: biological group assessment | `Figure_4_gene_group_biology.pdf` / `.png` | `redesign/complementarity_analysis/06_make_figures.R` | GO, Hallmark, Open Targets, and STRING group measurements compared with matched random sets |

Biological measurements are plotted as separate axes. STRING is secondary when
retained. The predictive ranking is calculated independently of these axes.

## Figure acceptance checks

The acceptance checks for the committed files are recorded in
`manifest.md`. They include the archived dataset set, finite numeric source
values, one-page PDF output, 300-dpi PNG output, and visual inspection for
clipping and label overlap. The rendered files and their SHA-256 hashes are
listed there.
