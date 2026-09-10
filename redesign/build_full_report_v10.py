#!/usr/bin/env python3
"""Build the corrected descriptive report for the benchmark datasets.

The repeated cross-validation splits share samples and training observations.
This report therefore gives fold-aggregated estimates and does not calculate
hypothesis-test p-values or binary win calls. All available methods are shown.
"""

import csv
import os
from collections import defaultdict
from pathlib import Path


CACHE_VERSION = "2026-08-26-corrected-v1"
PRIMARY_K = (10, 20, 50)
EXPECTED_SPLITS = 15
LABELS = {
    "GSE65682": "Sepsis, MARS, 28-day mortality",
    "GSE69683": "Asthma, U-BIOPRED, severe versus moderate",
    "GSE13355": "Psoriasis, paired positive control",
    "GSE107994": "Tuberculosis, active versus latent",
    "GSE101794": "Crohn's disease, RISK ileum",
    "imvigor210": "IMvigor210, anti-PD-L1 response",
    "sosall": "SOS-ALL, atopic dermatitis",
}


def result_root():
    return Path(os.environ.get(
        "GENESELECTR_REDESIGN_RESULTS_ROOT", "redesign/results_corrected"
    ))


def require_manifest(directory):
    path = directory / "run_manifest.version"
    if not path.exists():
        raise RuntimeError(f"missing corrected-run manifest: {directory}")
    version = path.read_text(encoding="utf-8").strip()
    if version != CACHE_VERSION:
        raise RuntimeError(f"incompatible run manifest: {directory}")


def require_extension_manifest(directory, evaluation_path):
    if evaluation_path.name == "eval_results.csv":
        return
    extension = evaluation_path.stem.removeprefix("eval_")
    path = directory / f"extension_manifest_{extension}.version"
    if not path.exists():
        raise RuntimeError(
            f"missing extension provenance for {evaluation_path.name}"
        )
    if path.read_text(encoding="utf-8").strip() != CACHE_VERSION:
        raise RuntimeError(
            f"incompatible extension provenance for {evaluation_path.name}"
        )


def normalise_arm(arm):
    if arm == "full_ungrouped":
        return "GS_full_ungrouped"
    if arm == "full_grouped":
        return "GS_full_grouped"
    return arm


def load_evaluations(directory):
    """Read all evaluation tables and reject conflicting result cells."""
    require_manifest(directory)
    values = {}
    fallbacks = {}
    sources = defaultdict(set)
    # The benchmark prediction-first table carries the same matched
    # Random rows as eval_results.csv. Read the base table first so an exact
    # repeated reference row retains the canonical source.
    paths = sorted(
        directory.glob("eval*.csv"),
        key=lambda path: (path.name != "eval_results.csv", path.name),
    )
    if not paths:
        raise RuntimeError(f"no evaluation tables found: {directory}")

    for path in paths:
        require_extension_manifest(directory, path)
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                required = {"repeat_idx", "fold_idx", "arm", "k", "AUC"}
                if not required <= row.keys():
                    continue
                arm = normalise_arm(row["arm"])
                key = (int(row["repeat_idx"]), int(row["fold_idx"]),
                       arm, int(row["k"]))
                value = float(row["AUC"])
                if key in values:
                    if values[key] != value:
                        raise RuntimeError(
                            f"conflicting evaluation cell for {key}: {path}"
                        )
                    sources[arm].add(path.name)
                    continue
                values[key] = value
                sources[arm].add(path.name)
                if "gate_fallback" in row and row["gate_fallback"] != "":
                    fallbacks[key] = row["gate_fallback"].lower() == "true"
    return values, fallbacks, sources


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def summaries(values, fallbacks):
    grouped = defaultdict(list)
    for (_, _, arm, k), auc in values.items():
        grouped[(arm, k)].append(auc)

    random = {k: mean(v) for (arm, k), v in grouped.items()
              if arm == "Random"}
    if not random:
        raise RuntimeError("matched Random results are missing")
    for cell, values_at_k in grouped.items():
        if len(values_at_k) != EXPECTED_SPLITS:
            raise RuntimeError(
                f"incomplete evaluation for {cell}: "
                f"{len(values_at_k)}/{EXPECTED_SPLITS} split rows"
            )
    methods = sorted({key[2] for key in values if key[2] != "Random"})
    rows = []
    for arm in methods:
        ks = sorted(k for (method, k) in grouped if method == arm)
        if ks != sorted(random):
            raise RuntimeError(
                f"panel sizes for {arm} do not match Random: {ks}"
            )
        deltas = []
        for k in PRIMARY_K:
            if (arm, k) not in grouped or k not in random:
                raise RuntimeError(
                    f"primary panel size k={k} is missing for {arm}"
                )
            deltas.append(mean(grouped[(arm, k)]) - random[k])
        arm_fallback = [flag for (_, _, method, _), flag in fallbacks.items()
                        if method == arm]
        rows.append({
            "arm": arm,
            "primary_delta": mean(deltas),
            "fallback_rate": mean([float(x) for x in arm_fallback])
                             if arm_fallback else 0.0,
            "by_k": [(k, mean(grouped[(arm, k)]),
                       mean(grouped[(arm, k)]) - random.get(k, float("nan")),
                       len(grouped[(arm, k)])) for k in ks],
        })
    return sorted(rows, key=lambda row: row["primary_delta"], reverse=True), random


def dataset_directories(root):
    directories = []
    validation = root / "validation_benchmark"
    if validation.exists():
        directories.extend((path.name, path) for path in sorted(validation.iterdir())
                           if path.is_dir())
    full_recipe = root / "full_recipe"
    if full_recipe.exists():
        directories.extend((path.name, path) for path in sorted(full_recipe.iterdir())
                           if path.is_dir())
    return directories


def render_dataset(name, directory):
    values, fallbacks, sources = load_evaluations(directory)
    rows, random = summaries(values, fallbacks)
    lines = [f"## {LABELS.get(name, name)}", "",
             f"Source directory: `{directory}`", "",
             "Primary delta is the arithmetic mean of AUC minus matched Random "
             "at k = 10, 20, and 50. Estimates are descriptive.", "",
             "| method | primary delta | gate fallback | evaluated files |",
             "|---|---:|---:|---|"]
    for row in rows:
        source_text = ", ".join(sorted(sources[row["arm"]]))
        lines.append(
            f"| {row['arm']} | {row['primary_delta']:+.3f} | "
            f"{100 * row['fallback_rate']:.1f}% | {source_text} |"
        )
    lines.extend(["", "### AUC by panel size", "",
                  "| method | k | mean AUC | AUC minus Random | split rows |",
                  "|---|---:|---:|---:|---:|"])
    for row in rows:
        for k, auc, delta, count in row["by_k"]:
            lines.append(
                f"| {row['arm']} | {k} | {auc:.3f} | {delta:+.3f} | {count} |"
            )
    lines.extend(["", "Matched Random mean AUC: " + ", ".join(
        f"k={k}: {auc:.3f}" for k, auc in sorted(random.items())), ""])
    return lines


def main():
    root = result_root()
    directories = dataset_directories(root)
    if not directories:
        raise RuntimeError(f"no corrected benchmark directories found under {root}")

    lines = ["# GeneSelectR benchmark report", "",
             "## Summary", "",
             "This report includes every method found in the corrected evaluation "
             "tables. SOS-ALL was used for method development. The other six "
             "datasets were used for validation after development.", "",
             "Repeated cross-validation folds and repeats share observations. The "
             "tables report descriptive means. No split-level hypothesis tests or "
             "equivalence claims are calculated.", ""]
    for name, directory in directories:
        lines.extend(render_dataset(name, directory))

    output = root / "full_benchmark" / "REPORT.md"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
