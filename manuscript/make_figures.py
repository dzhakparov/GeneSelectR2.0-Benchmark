#!/usr/bin/env python3
"""Render manuscript Figures 2-4 and the biological-enrichment counts from the saved source tables."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42})

LABEL = {"GS_full_ungrouped": "GeneSelectR 2.0", "Stabl": "Stabl",
         "mRMR": "mRMR", "RF_importance": "Random forest",
         "DGE": "DGE", "Boruta": "Boruta", "LASSO": "LASSO",
         "ElasticNet": "Elastic net"}
COL = {"GS_full_ungrouped": "#006C80", "Stabl": "#A64279",
       "mRMR": "#52835E", "RF_importance": "#C06739", "DGE": "#777777",
       "Boruta": "#8060A5", "LASSO": "#B38721", "ElasticNet": "#4678B7"}
DATASETS = ["GSE101794", "GSE107994", "GSE13355", "GSE65682",
            "GSE69683", "imvigor210", "sosall"]
DS_LABEL = ["Crohn's disease", "Tuberculosis", "Psoriasis", "Sepsis",
            "Asthma", "IMvigor210", "SOS-ALL"]

def read(name):
    return pd.read_csv(DATA / name)

def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

# Figure 1 is drawn by make_workflow_figure.py.

# Figure 2: primary predictive comparison, including Stabl.
primary = read("by_dataset_primary.csv")
validation = primary[primary.dataset != "sosall"].groupby("method").AUC.mean().sort_values(ascending=False)
development = primary[primary.dataset == "sosall"].set_index("method").AUC
assert len(validation) == len(development) == 8
order = validation.index.tolist()
fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8), layout="constrained", sharex=True)
for ax, vals, title in zip(axes, [validation, development],
                           ["A  Six validation datasets", "B  SOS-ALL development dataset"]):
    for i, m in enumerate(order):
        ax.scatter(vals[m], i, color=COL[m], marker="o", s=42, zorder=3)
        ax.text(vals[m]+.006, i, f"{vals[m]:.4f}", va="center", fontsize=7.8)
    ax.set(yticks=range(8), yticklabels=[LABEL[m] for m in order], xlim=(.58,.96),
           xlabel="Mean outer-test AUC")
    ax.invert_yaxis(); ax.grid(axis="x", color="#E5E9EE")
    ax.set_title(title, loc="left", weight="bold", fontsize=10)
save(fig,"Figure_2_prediction")

# Figure 3: Nogueira stability across the 15 outer selected sets.
old = read("stability_classical.csv")
old = old[old.dataset.isin(DATASETS) & old.k.isin([10,20,50])]
old.method = old.method.replace({"GeneSelectR":"GS_full_ungrouped",
                                 "Random forest":"RF_importance","Elastic net":"ElasticNet"})
old = old.groupby(["dataset","method"],as_index=False).nogueira_stability.mean()
new = read("stability_by_dataset_primary.csv")
stabl = new[new.method == "Stabl"].rename(columns={"Nogueira":"nogueira_stability"})
stab = pd.concat([old,stabl],ignore_index=True)
assert len(stab) == 7*8 and not stab.duplicated(["dataset","method"]).any()
fig, axes = plt.subplots(2,4,figsize=(10,5.5),layout="constrained",sharex=True)
for ax,ds,title in zip(axes.flat,DATASETS,DS_LABEL):
    d=stab[stab.dataset==ds].set_index("method").loc[order]
    for i,m in enumerate(order):
        ax.scatter(d.loc[m,"nogueira_stability"],i,color=COL[m],marker="o",s=24,zorder=3)
        ax.annotate(f"{d.loc[m,'nogueira_stability']:.3f}",
                    xy=(d.loc[m,"nogueira_stability"],i),
                    xytext=(6,0),textcoords="offset points",
                    ha="left",va="center",fontsize=7)
    ax.set(yticks=range(8),yticklabels=[LABEL[m] for m in order],xlim=(0,.95),
           xticks=[0,.25,.5,.75])
    ax.invert_yaxis();ax.grid(axis="x",color="#E5E9EE")
    ax.set_title(title,loc="left",weight="bold",fontsize=10)
axes[1,3].axis("off")
fig.supxlabel("Nogueira selection stability across 15 outer divisions",fontsize=10)
save(fig,"Figure_3_selection_stability")

# Figure 4: asthma gene-level evidence plus k=20 gene-set biological assessment.
eight = read("asthma_eight.csv").sort_values(["n_outer_top20","gene"],ascending=[False,True])
asthma_auc = read("by_dataset_size.csv").query("dataset == 'GSE69683' and k == 20").set_index("method").AUC.sort_values(ascending=False)
paired=read("asthma_k20_biology_comparison.csv")
bio_labels=["GO-BP","Hallmark","STRING","OT 0.05","OT 0.10"]
plot_methods=[LABEL[m] for m in asthma_auc.index]
assert len(eight)==8 and len(asthma_auc)==8 and len(paired)==8*5
assert set(paired.source)==set(bio_labels) and set(paired.method)==set(plot_methods)
bio_matrix=(paired.pivot(index="method",columns="source",values="median_enrichment_ratio")
                  .loc[plot_methods,bio_labels])
assert bio_matrix.notna().all().all()
gs_summary=paired[paired.method=="GeneSelectR 2.0"].set_index("source").loc[bio_labels].reset_index()
gs_summary[["source","median_enrichment_ratio","outer_divisions"]].to_csv(
    DATA/"asthma_k20_biology_summary.csv",index=False)
y=np.arange(8);names=eight.gene.tolist()
fig=plt.figure(figsize=(9.6,10.3),layout="constrained")
grid=fig.add_gridspec(3,2,height_ratios=[1.15,1.05,1.45])
ax=fig.add_subplot(grid[0,0])
for i,(m,v) in enumerate(asthma_auc.items()):
    ax.scatter(v,i,color=COL[m],marker="o",s=35)
    ax.text(v+.0008,i,f"{v:.3f}",va="center",fontsize=7.5)
ax.set(yticks=y,yticklabels=[LABEL[m] for m in asthma_auc.index],xlim=(.765,.807),
       xlabel="Mean outer-test AUC")
ax.invert_yaxis();ax.grid(axis="x",color="#E5E9EE")
ax.set_title("A  Prediction with 20 genes",loc="left",weight="bold",fontsize=10)
ax=fig.add_subplot(grid[0,1])
ax.barh(y-.17,eight.n_outer_top20,height=.32,color=COL["GS_full_ungrouped"],label="GeneSelectR 2.0")
ax.barh(y+.17,eight.DGE_top20_count,height=.32,color="#A0A5AA",label="DGE")
ax.set(yticks=y,yticklabels=names,xlim=(0,16),xticks=[0,5,10,15],xlabel="Outer selections (of 15)")
ax.invert_yaxis()
ax.set_title("B  Repeated gene selection",loc="left",weight="bold",fontsize=10)
ax.set_ylim(7.6,-1.5)
ax.legend(loc="upper right",ncol=2,frameon=False,fontsize=8,handlelength=1.2,borderaxespad=0.1)
ax=fig.add_subplot(grid[1,0])
ax.scatter(eight.mean_geneselectr_calibrated_utility,y,color=COL["GS_full_ungrouped"],s=32)
ax.axvline(1,color="#88909A",ls="--",lw=.8)
ax.set(yticks=y,yticklabels=names,xlim=(0,16),xlabel="Mean adjusted utility ratio")
ax.invert_yaxis();ax.set_title("C  Predictive utility",loc="left",weight="bold",fontsize=10)
ax=fig.add_subplot(grid[1,1])
for i,(_,row) in enumerate(eight.iterrows()):
    if pd.notna(row.OT_association_score):
        ax.scatter(row.OT_association_score,i,color="#8060A5",s=32)
        ax.text(row.OT_association_score+.006,i,f"{row.OT_association_score:.4f}",va="center",fontsize=7)
    else:
        ax.text(.008,i,"No recorded score",va="center",fontsize=7,color="#68727E")
ax.set(yticks=y,yticklabels=names,xlim=(0,.43),xlabel="Open Targets asthma association score")
ax.invert_yaxis();ax.set_title("D  Gene-level disease association",loc="left",weight="bold",fontsize=10)
ax=fig.add_subplot(grid[2,:])
log_ratio=np.log2(np.clip(bio_matrix.to_numpy(),1/16,16))
im=ax.pcolormesh(log_ratio,cmap="RdBu_r",
                 norm=TwoSlopeNorm(vmin=-4,vcenter=0,vmax=4),
                 edgecolors="white",linewidth=.35)
ax.set(xticks=np.arange(len(bio_labels))+.5,xticklabels=bio_labels,
       yticks=np.arange(len(plot_methods))+.5,yticklabels=plot_methods,
       xlim=(0,len(bio_labels)),ylim=(len(plot_methods),0))
for i in range(len(plot_methods)):
    for j in range(len(bio_labels)):
        value=bio_matrix.iloc[i,j]
        colour="white" if abs(log_ratio[i,j])>=2.25 else "black"
        ax.text(j+.5,i+.5,f"{value:.2f}",ha="center",va="center",fontsize=7.3,color=colour)
ax.set_title("E  Biology of the selected 20-gene sets",loc="left",weight="bold",fontsize=10)
colourbar=fig.colorbar(im,ax=ax,fraction=.025,pad=.025)
colourbar.set_label("Log2 enrichment ratio",fontsize=8)
colourbar.ax.tick_params(labelsize=7)
colourbar.solids.set_rasterized(False)
save(fig,"Figure_4_asthma_case")

# Supplementary figure: cross-dataset biological enrichment counts for all eight methods (Section 4.3).
comparison = read("biology_comparison_by_dataset.csv")
assert comparison.shape == (7 * 8, 7)
assert not comparison.duplicated(["dataset", "method"]).any()
method_labels = [LABEL[m] for m in order]
assert set(comparison.method) == set(method_labels)
assert set(comparison.dataset) == set(DATASETS)
enrichment_counts = (
    comparison.melt(
        id_vars=["dataset", "method"], value_vars=bio_labels,
        var_name="source", value_name="enrichment_ratio"
    )
    .assign(above_reference=lambda frame: frame.enrichment_ratio > 1)
    .groupby(["method", "source"], as_index=False, observed=False)
    .above_reference.sum()
)
count_matrix = (
    enrichment_counts.pivot(index="method", columns="source", values="above_reference")
    .loc[method_labels, bio_labels]
    .astype(int)
)
assert count_matrix.shape == (8, 5)
enrichment_counts = (
    count_matrix.rename_axis(index="method", columns="source")
    .stack()
    .rename("datasets_above_reference")
    .reset_index()
)
enrichment_counts.to_csv(DATA / "biology_enrichment_counts.csv", index=False)

label_to_code = {label: code for code, label in LABEL.items()}
fig, axes = plt.subplots(
    1, 5, figsize=(10.2, 3.4), layout="constrained", sharex=True, sharey=True
)
for panel, (ax, source) in enumerate(zip(axes, bio_labels)):
    for i, method in enumerate(method_labels):
        value = int(count_matrix.loc[method, source])
        ax.scatter(value, i, color=COL[label_to_code[method]], marker="o", s=34, zorder=3)
        ax.text(value + .18, i, str(value), va="center", fontsize=7.2)
    ax.set(
        xlim=(-.2, 7.75), xticks=[0, 2, 4, 6],
        yticks=np.arange(len(method_labels)), ylim=(len(method_labels)-.5, -.5),
    )
    ax.grid(axis="x", color="#E5E9EE")
    ax.set_title(f"{chr(65 + panel)}  {source}", loc="left", weight="bold", fontsize=9)
    if panel == 0:
        ax.set_yticklabels(method_labels)
    else:
        ax.tick_params(axis="y", left=False, labelleft=False)
fig.supxlabel("Datasets with enrichment ratio above 1 (of 7)", fontsize=9)
save(fig, "Supplementary_biological_enrichment_counts")

# A single table joins the four source-specific dataset-level summaries.
biology=read("biology_by_dataset.csv").query("method == 'GS_full_ungrouped'").copy()
ots=read("ot_by_dataset.csv").query("method == 'GS_full_ungrouped'")
wide=ots.pivot(index="dataset",columns="cutoff",values="open_targets_enrichment")
wide.columns=[f"open_targets_{v:.2f}" for v in wide.columns]
table=biology.set_index("dataset")[["GO_semantic_enrichment","hallmark_enrichment","string_enrichment"]].join(wide)
table=table.loc[DATASETS].reset_index()
assert table.shape==(7,6) and table.notna().all().all()
table.to_csv(DATA/"biology_dataset_summary.csv",index=False)
print("Rendered Figures 2-4 and the supplementary enrichment-count figure; all source-data checks passed.")
