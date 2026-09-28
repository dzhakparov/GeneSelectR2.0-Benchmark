#!/usr/bin/env python3
"""Vector GeneSelectR 2.0 workflow figure (Figure 1). Writes PDF, SVG and PNG."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle
from pathlib import Path

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "path"})
NAVY = "#1F2A44"; GREY = "#6B7280"
BLUE, BLUE_E = "#EEF4FB", "#6A9FD8"
GREEN, GREEN_E = "#EEF6EA", "#8CC47A"
PEACH, PEACH_E = "#FDF0E8", "#EE8F63"
CASE, CTRL = "#D9687F", "#4E8FD6"

fig = plt.figure(figsize=(12, 7.3))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 120); ax.set_ylim(0, 73); ax.axis("off")

def panel(x, y, w, h, fc, ec, title):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.6",
                                fc=fc, ec=ec, lw=1.4, zorder=1))
    ax.text(x + w / 2, y + h - 3.2, title, ha="center", va="center", fontsize=12.5,
            weight="bold", color=NAVY, zorder=5)

def footer(x, y, w, lines):
    ax.text(x + w / 2, y + 4.2, "\n".join(lines), ha="center", va="center", fontsize=9.6,
            color=NAVY, linespacing=1.45, zorder=5)

def person(cx, cy, s, c):
    ax.add_patch(Circle((cx, cy + 1.05 * s), 0.42 * s, fc=c, ec="none", zorder=4))
    ax.add_patch(FancyBboxPatch((cx - 0.55 * s, cy - 0.9 * s), 1.1 * s, 1.35 * s,
                                boxstyle=f"round,pad=0,rounding_size={0.35*s}", fc=c, ec="none", zorder=4))

def arrow(x0, y0, x1, y1, lw=2.4, c="#8C939D", ms=18):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms,
                                 lw=lw, color=c, zorder=3))

def inset(x, y, w, h):
    """Axes placed in data coordinates of the main canvas."""
    return fig.add_axes([x / 120, y / 73, w / 120, h / 73])

# geometry
PW, PH = 36.5, 29.5
XS = [2, 41.75, 81.5]; YT, YB = 39.0, 2.5

# bracket "Training data only"
ax.plot([XS[0] + 0.5, XS[0] + 0.5, 51.5], [70.2, 71.2, 71.2], color=BLUE_E, lw=1.3)
ax.plot([68.5, XS[2] + PW - 0.5, XS[2] + PW - 0.5], [71.2, 71.2, 70.2], color=BLUE_E, lw=1.3)
ax.text(60, 71.2, "Training data only", ha="center", va="center", fontsize=11, color="#2F6DB5")

# ---------- Panel 1: training expression and outcome
x, y = XS[0], YT
panel(x, y, PW, PH, BLUE, BLUE_E, "Training expression and outcome")
rng = np.random.default_rng(3)
hm = rng.normal(0, 0.55, (8, 10)); hm[:, 3:7] += np.array([0.9, 1.5, 1.7, 1.0])
a = inset(x + 3.4, y + 10.2, 17.5, 13.2)
a.imshow(hm, cmap="RdBu_r", vmin=-2.3, vmax=2.3, aspect="auto")
a.set_xticks([]); a.set_yticks([])
for s in a.spines.values(): s.set_color("#9AA3AD")
for k in range(1, 10): a.axvline(k - 0.5, color="white", lw=0.8)
for k in range(1, 8): a.axhline(k - 0.5, color="white", lw=0.8)
a.set_ylabel("Samples", fontsize=9.2, color=NAVY, labelpad=3); a.set_xlabel("Genes", fontsize=9.2, color=NAVY, labelpad=3)
ax.plot([x + 23.6, x + 23.6], [y + 9, y + 24.5], color="#C6D3E3", lw=1)
ax.text(x + 30, y + 23.2, "Outcome", ha="center", fontsize=9.6, weight="bold", color=NAVY)
for i in range(4):
    person(x + 25.9 + 2.75 * i, y + 18.6, 1.25, CTRL)
    person(x + 25.9 + 2.75 * i, y + 14.0, 1.25, CASE)
ax.add_patch(Circle((x + 26.2, y + 10.9), 0.6, fc=CASE, ec="none")); ax.text(x + 27.3, y + 10.9, "Case", va="center", fontsize=9, color=NAVY)
ax.add_patch(Circle((x + 26.2, y + 8.7), 0.6, fc=CTRL, ec="none")); ax.text(x + 27.3, y + 8.7, "Control", va="center", fontsize=9, color=NAVY)
footer(x, y, PW, ["Train-only processing", "2,000 variable genes"])

# ---------- Panel 2: repeated elastic-net models
x = XS[1]
panel(x, y, PW, PH, BLUE, BLUE_E, "Repeated elastic-net models")
ax.text(x + 8.6, y + 22.0, "10 repeats\nof 5 folds", ha="center", va="center", fontsize=9.4, color=NAVY, linespacing=1.3)
rows = [y + 17.6, y + 15.2, y + 12.8, y + 9.6]
for r, yy in enumerate(rows):
    for c in range(5):
        hl = (c == [0, 1, 2, 4][r])
        ax.add_patch(Rectangle((x + 2.6 + c * 2.45, yy), 2.25, 1.55, fc="#6FA8E6" if hl else "#C9DDF3", ec="none", zorder=3))
ax.text(x + 8.6, y + 11.5, "⋮", ha="center", va="center", fontsize=11, color=GREY)
ax.plot([x + 15.3, x + 15.9, x + 15.9, x + 15.3], [y + 19.3, y + 19.3, y + 9.4, y + 9.4], color=NAVY, lw=1)
arrow(x + 15.9, y + 14.4, x + 18.3, y + 14.4, lw=1.1, c=NAVY, ms=9)
a = inset(x + 21.3, y + 9.3, 13.2, 14.8)
lam = np.linspace(0, 1, 100)
for k, (c0, col) in enumerate(zip([0.95, 0.7, 0.45, 0.25, -0.2, -0.45, -0.75],
                                    ["#4E8FD6", "#D9687F", "#E0453A", "#7FB069", "#9D7CC6", "#57A773", "#8C939D"])):
    a.plot(lam, c0 * np.clip(1 - lam / (0.55 + 0.06 * k), 0, None) ** 1.3, color=col, lw=1.3)
a.axvline(0.62, color=GREY, ls="--", lw=0.9)
a.set_xticks([]); a.set_yticks([]); a.spines[["top", "right"]].set_visible(False)
a.set_xlabel(r"log($\lambda$)", fontsize=8.6, color=NAVY, labelpad=2); a.set_ylabel("Coefficients", fontsize=8.6, color=NAVY, labelpad=2)
a.set_title("Elastic net", fontsize=9, color=NAVY, pad=3); a.patch.set_alpha(0)
footer(x, y, PW, ["50 fits: 10 repeats of 5 folds", "Recurrence and excluded-sample utility"])

# ---------- Panel 3: outcome-permutation comparison
x = XS[2]
panel(x, y, PW, PH, BLUE, BLUE_E, "Outcome-permutation comparison")
ax.text(x + 8.4, y + 22.3, "Permute\noutcome", ha="center", va="center", fontsize=9.6, weight="bold", color=NAVY, linespacing=1.25)
cols_top = [CTRL, CASE, CTRL, CASE, CTRL]; cols_bot = [CASE, CTRL, CTRL, CASE, CASE]
for i in range(5):
    person(x + 3.4 + 2.5 * i, y + 17.0, 1.05, cols_top[i]); person(x + 3.4 + 2.5 * i, y + 9.7, 1.05, cols_bot[i])
arrow(x + 5.6, y + 15.3, x + 11.2, y + 12.7, lw=1.1, c=NAVY, ms=8); arrow(x + 11.2, y + 15.3, x + 5.6, y + 12.7, lw=1.1, c=NAVY, ms=8)
ax.plot([x + 16.4, x + 16.4], [y + 9, y + 24.5], color="#C6D3E3", lw=1)
a = inset(x + 19.2, y + 9.3, 15.6, 14.8)
g = np.linspace(-3, 7, 300)
a.fill_between(g, np.exp(-g ** 2 / 2), color="#B9C0C9", alpha=.55); a.plot(g, np.exp(-g ** 2 / 2), color=GREY, lw=1.3, label="Permuted")
a.plot(g, np.exp(-(g - 3.8) ** 2 / 2), color="#D23F4C", lw=1.4, label="Observed")
a.axvline(0, color=GREY, ls=":", lw=0.8)
a.set_xticks([]); a.set_yticks([]); a.spines[["top", "right"]].set_visible(False); a.set_ylim(0, 1.35)
a.set_xlabel("Component score", fontsize=8.6, color=NAVY, labelpad=2); a.set_ylabel("Density", fontsize=8.6, color=NAVY, labelpad=2)
a.legend(frameon=False, fontsize=7.8, loc="upper right", handlelength=1.3, borderaxespad=0.1); a.patch.set_alpha(0)
footer(x, y, PW, ["20 permutations, 20 fits each", "One ratio for each component"])

# ---------- Panel 6: rank genes (bottom right)
x, y = XS[2], YB
panel(x, y, PW, PH, GREEN, GREEN_E, "Rank genes")
labels = ["Gene 1", "Gene 2", "Gene 3", "Gene 4", "⋮", "Gene N"]; vals = [1, .82, .68, .55, None, .12]
for i, (lab, v) in enumerate(zip(labels, vals)):
    yy = y + 22.3 - i * 2.35
    ax.text(x + 7.6, yy, lab, ha="right", va="center", fontsize=8.8, color=NAVY)
    if v is not None:
        ax.add_patch(Rectangle((x + 8.3, yy - 0.75), 11.5 * v, 1.5, fc="#7FBF6A" if i < 4 else "#BFE0B4", ec="none", zorder=3))
ax.plot([x + 8.3, x + 8.3], [y + 9.0, y + 23.6], color=NAVY, lw=1)
ax.text(x + 28.5, y + 18.6, "Rank score", ha="center", va="center", fontsize=10.5, color=NAVY)
ax.text(x + 28.5, y + 14.7, r"$\sqrt{r_1 \times r_2}$", ha="center", va="center", fontsize=15, color=NAVY)
footer(x, y, PW, ["Geometric mean of the two ratios", "Select a fixed or training-chosen size"])

# ---------- Panel 5: evaluate the selected gene set
x = XS[1]
panel(x, y, PW, PH, GREEN, GREEN_E, "Evaluate the selected gene set")
ax.text(x + 7.4, y + 22.3, "Outer test\nsamples", ha="center", va="center", fontsize=9.6, weight="bold", color=NAVY, linespacing=1.25)
for i in range(4): person(x + 3.6 + 2.6 * i, y + 16.0, 1.2, "#8C939D")
ax.text(x + 7.4, y + 10.9, "Held-out\nsamples", ha="center", va="center", fontsize=9, color=GREY, linespacing=1.25)
ax.plot([x + 14.6, x + 14.6], [y + 9, y + 24.5], color="#CFE3C7", lw=1)
a = inset(x + 19.6, y + 10.6, 14.3, 13.6)
fpr = np.linspace(0, 1, 200); tpr = 1 - (1 - fpr) ** 3.1
a.plot(fpr, tpr, color="#2F6DB5", lw=1.6); a.plot([0, 1], [0, 1], color=GREY, ls="--", lw=0.8)
a.set_xlim(0, 1); a.set_ylim(0, 1); a.set_xticks([0, .5, 1]); a.set_yticks([0, .5, 1])
a.tick_params(labelsize=7.2, length=2, colors=NAVY)
a.set_xlabel("1 − Specificity", fontsize=8.6, color=NAVY, labelpad=1); a.set_ylabel("Sensitivity", fontsize=8.6, color=NAVY, labelpad=1)
a.text(0.95, 0.1, "AUC = 0.78", ha="right", fontsize=8.4, color=NAVY); a.patch.set_alpha(0)
footer(x, y, PW, ["Common prediction procedure", "Held-out outer samples"])

# ---------- Panel 4: interpret selected genes
x = XS[0]
panel(x, y, PW, PH, PEACH, PEACH_E, "Interpret selected genes")
pts = np.array([[0, .5], [.3, .9], [.35, .2], [.65, .6], [.95, .9], [1, .25], [.62, .05]])
P = np.c_[x + 3.2 + pts[:, 0] * 11.5, y + 12.2 + pts[:, 1] * 10.2]
for i, j in [(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (3, 5), (2, 6), (5, 6), (1, 4)]:
    ax.plot(*P[[i, j]].T, color="#8FA9C9", lw=1.6, zorder=3)
for k, (px, py) in enumerate(P):
    ax.add_patch(Circle((px, py), 0.95, fc=["#3D7FD0", "#3D7FD0", "#6FA3E0", CASE, "#6FA3E0", "#3D7FD0", "#6FA3E0"][k], ec="white", lw=1.2, zorder=4))
ax.text(x + 9.0, y + 9.6, "Gene-set annotation", ha="center", va="center", fontsize=9, color=NAVY)
ax.plot([x + 18.2, x + 18.2], [y + 9, y + 24.5], color="#F2C8B3", lw=1)
a = inset(x + 20.2, y + 12.2, 14.6, 10.6)
g = np.linspace(-3, 6.4, 300); a.fill_between(g, np.exp(-g ** 2 / 2), color="#B9C0C9", alpha=.9)
a.plot(g, np.exp(-g ** 2 / 2), color=GREY, lw=1.1); a.plot([4.3, 4.3], [0, 0.8], color="#C94F6D", lw=2)
a.plot([4.3], [0.8], "o", color="#C94F6D", ms=5); a.axhline(0, color=GREY, lw=1)
a.axis("off"); a.set_ylim(-0.02, 1.05)
ax.text(x + 27.5, y + 9.6, "Enrichment vs random", ha="center", va="center", fontsize=9, color=NAVY)
footer(x, y, PW, ["GO, Hallmark, Open Targets", "Compare with matched random sets"])

# ---------- badges
def badge(cx, cy, text, fc, ec, tc):
    ax.text(cx, cy, text, ha="center", va="center", fontsize=9.2, color=tc, linespacing=1.25, zorder=6,
            bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.9", fc=fc, ec=ec, lw=1))
badge(XS[0] + PW / 2, 34.7, "Biological annotation\ndoes not alter ranking", "#FBE1D6", "#F4B79C", "#B03A2E")
badge(XS[1] + PW / 2, 34.7, "Outer test samples quantify\npredictive performance", "#DDEFD6", "#A9D59A", "#2E6B2A")

# ---------- flow arrows
arrow(XS[0] + PW + 0.3, YT + PH / 2, XS[1] - 0.3, YT + PH / 2)
arrow(XS[1] + PW + 0.3, YT + PH / 2, XS[2] - 0.3, YT + PH / 2)
arrow(XS[2] + PW / 2, YT - 0.3, XS[2] + PW / 2, YB + PH + 0.3)
arrow(XS[2] - 0.3, YB + PH / 2, XS[1] + PW + 0.3, YB + PH / 2)
arrow(XS[1] - 0.3, YB + PH / 2, XS[0] + PW + 0.3, YB + PH / 2)

out = Path(__file__).resolve().parent / "figures"; out.mkdir(exist_ok=True)
fig.savefig(out / "Figure_1_workflow.svg")  # text stored as outlines (svg.fonttype = "path")
fig.savefig(out / "Figure_1_workflow.png", dpi=300)
# Vector PDF with outlined text, rendered from the SVG (requires cairosvg).
import cairosvg
cairosvg.svg2pdf(url=str(out / "Figure_1_workflow.svg"), write_to=str(out / "Figure_1_workflow.pdf"))
print("written", out)
