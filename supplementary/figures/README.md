# Supplementary figures

## Verification state

The clean archive worktree contains no ignored benchmark result tables or
figure source data. No rendered figure is committed in this branch. Figure
generation remains available through
`redesign/complementarity_analysis/06_make_figures.R` after scripts 01 to 04
have produced verified CSV inputs.

## Figure inventory

| Figure | Generation script | Source analysis | Caption |
|---|---|---|---|
| Figure 2: prediction and stability | `redesign/complementarity_analysis/06_make_figures.R` | Saved prediction and gene-set stability tables | Mean test AUC and cross-split gene-set stability by method, dataset group, and panel size |
| Supplementary stability curves and heatmap | `redesign/complementarity_analysis/06_make_figures.R` | `gene_set_stability_summary.csv` and recurrence tables | Stability measures across panel sizes and datasets |
| Figure 3: DGE/GeneSelectR complementarity | `redesign/complementarity_analysis/06_make_figures.R` | DGE/GeneSelectR concordance and recurrence tables | Rank percentiles, top-k overlap, and recurrent genes from matched candidate pools |
| Supplementary DGE score components | `redesign/complementarity_analysis/06_make_figures.R` | Recomputed DGE and GeneSelectR component scores | DGE rank compared with GeneSelectR score components |
| Figure 4: biological group assessment | `redesign/complementarity_analysis/06_make_figures.R` | Gene-group biology tables | GO, Hallmark, Open Targets, and STRING group measurements compared with matched random sets |

Biological measurements are plotted as separate axes. STRING is secondary when
retained. The predictive ranking is calculated independently of these axes.

## Figure acceptance checks

Before a rendered figure is accepted for the supplement, verify the source CSV,
dataset role, gene and gene-set terminology, panel-size labels, and output
dimensions. Confirm that all figure outputs have a corresponding source table
and generation script. Store accepted PNG, SVG, or PDF files directly in this
directory after the checks pass.
