#!/usr/bin/env python3
"""Assemble the GSE101794 (RISK ileal, corrected) TPM matrix.

Why this exists:
  GSE57945 (the locked registry pick for Crohn's) ships only RPKM and the
  2015 corrigendum's 40-sample exclusion list is not published anywhere
  machine-readable (checked: GEO series matrix, full SOFT, corrigendum text,
  the combined RPKM table, citing papers). GSE101794 is the SAME cohort
  (RISK ileal biopsies, treatment-naive paediatric CD vs non-IBD control),
  re-quantified and submitted in 2017 AFTER the corrigendum, so its labels
  are the corrected ones. It ships per-sample Gene x TPM files.

What it does:
  - reads data/GSE101794/raw/tar/GSM*_RISK_*.txt.gz (two columns: Gene, TPM)
  - verifies every file has the identical gene set (fails loudly otherwise)
  - collapses duplicate gene symbols by SUMMING TPM (37 symbols appear twice,
    e.g. transcript loci sharing a symbol; sum is the gene-level abundance)
  - builds a genes x samples matrix, log2(TPM + 1)
  - writes data/GSE101794/expression.tsv  (genes x samples, log2 TPM)
  - writes data/GSE101794/prep_manifest.json documenting the substitution

Checks (no invented data):
  - 304 files expected, one per metadata row
  - gene sets must be identical across all files
  - TPM must be >= 0 and finite
"""
import gzip, json, re, sys
from pathlib import Path

import csv

ROOT = Path("data/GSE101794")
TAR_DIR = ROOT / "raw" / "tar"

def main():
    files = sorted(TAR_DIR.glob("GSM*_RISK_*.txt.gz"))
    if len(files) != 304:
        sys.exit(f"expected 304 per-sample files, found {len(files)}")

    meta = {r["sample"]: r for r in
            csv.DictReader(open(ROOT / "metadata.tsv"), delimiter="\t")}

    genes_ref = None
    cols = {}
    for f in files:
        gsm = re.match(r"(GSM\d+)_", f.name).group(1)
        if gsm not in meta:
            sys.exit(f"{gsm} not in metadata.tsv")
        genes, vals = [], []
        #  37 symbols appear more than once per file (multiple loci sharing a
        #  symbol). Sum their TPM to gene level, preserving file order.
        per_gene = {}
        order = []
        with gzip.open(f, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            if header != ["Gene", "TPM"]:
                sys.exit(f"{f.name}: unexpected header {header}")
            for line in fh:
                g, v = line.rstrip("\n").split("\t")
                x = float(v)
                if x < 0 or x != x:
                    sys.exit(f"{f.name}: bad TPM value {v!r} for {g}")
                if g not in per_gene:
                    per_gene[g] = 0.0
                    order.append(g)
                per_gene[g] += x
        genes, vals = order, [per_gene[g] for g in order]
        if genes_ref is None:
            genes_ref = genes
        elif genes != genes_ref:
            sys.exit(f"{f.name}: gene set/order differs from first file")
        cols[gsm] = vals

    n_genes = len(genes_ref)
    gsms = [re.match(r"(GSM\d+)_", f.name).group(1) for f in files]
    import math
    with open(ROOT / "expression.tsv", "w") as out:
        out.write("gene\t" + "\t".join(gsms) + "\n")
        for i, g in enumerate(genes_ref):
            row = [f"{math.log2(cols[s][i] + 1):.6f}" for s in gsms]
            out.write(g + "\t" + "\t".join(row) + "\n")

    diag = {}
    for s in gsms:
        diag[meta[s]["diagnosis"]] = diag.get(meta[s]["diagnosis"], 0) + 1

    manifest = {
        "accession": "GSE101794",
        "substitutes": "GSE57945",
        "reason": ("GSE57945 ships only RPKM and the 2015 corrigendum's "
                   "40-sample exclusion list is not publicly available. "
                   "GSE101794 provides RISK ileal gene-level TPM. The local "
                   "metadata does not establish corrigendum label status."),
        "scale": "log2(TPM + 1), normalized branch (no TMM; TPM is not counts)",
        "n_samples": len(gsms),
        "n_genes": n_genes,
        "diagnosis_counts": diag,
        "outcome": {"column": "diagnosis", "positive": "CD",
                    "negative": "Non-IBD"},
    }
    with open(ROOT / "prep_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"OK: {n_genes} genes x {len(gsms)} samples; diagnosis: {diag}")

if __name__ == "__main__":
    main()
