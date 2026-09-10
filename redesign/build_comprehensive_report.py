#!/usr/bin/env python3
"""Build comprehensive predictive and biological redesign comparisons.

The benchmark uses three repeated five-fold cross-validation runs. Folds share
training observations across repeats. All uncertainty summaries are descriptive;
split-level hypothesis tests are omitted because their independence assumption
does not hold. Corrected-run and extension manifests are required in strict mode.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import os
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path


CACHE_VERSION = "2026-08-26-corrected-v1"
PRIMARY_K = (10, 20, 50)
PANEL_SIZES = (10, 20, 50, 100, 200, 500)
EXPECTED_SPLITS = 15
EXPECTED_METHODS = 27
EXCLUDED_METHODS = {"dc_pf_ctrl"}
DETERMINISTIC_EVALUATION_VERSION = "2026-08-29-deterministic-v1"
EXPECTED_EVALUATIONS = {
    "eval_results.csv",
    "eval_predfirst.csv",
    "eval_ensemble.csv",
    "eval_softprior.csv",
    "eval_datadriven.csv",
    "eval_horseshoe.csv",
    "eval_ensemble2.csv",
    "eval_wsweep.csv",
    "eval_hsstab.csv",
    "eval_prune.csv",
    "eval_adapt.csv",
}
LABELS = {
    "GSE65682": "Sepsis, MARS, 28-day mortality",
    "GSE69683": "Asthma, U-BIOPRED, severe versus moderate",
    "GSE13355": "Psoriasis, paired positive control",
    "GSE107994": "Tuberculosis, active versus latent",
    "GSE101794": "Crohn's disease, RISK ileum",
    "imvigor210": "IMvigor210, anti-PD-L1 response",
    "sosall": "SOS-ALL, atopic dermatitis",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        default=os.environ.get(
            "GENESELECTR_REDESIGN_RESULTS_ROOT", "redesign/results_corrected"
        ),
    )
    parser.add_argument("--output-dir")
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument(
        "--evaluation-file",
        help="Load only this evaluation CSV from each dataset directory.",
    )
    return parser.parse_args()


def normalise_arm(arm):
    return {
        "full_ungrouped": "GS_full_ungrouped",
        "full_grouped": "GS_full_grouped",
    }.get(arm, arm)


def method_family(arm):
    if arm in {"DGE", "LASSO", "ElasticNet", "mRMR", "Boruta", "RF_importance"}:
        return "Classical comparator"
    if arm.startswith("GS_full"):
        return "Base GeneSelectR"
    if arm.startswith("predfirst"):
        return "Prediction-first"
    if arm in {"soft_gs", "soft_pf", "cb_gs_w2"}:
        return "Biology-weighted"
    if arm in {"dc_gs", "dc_pf", "cb_gs", "cb_pf"}:
        return "Data-driven or combined prior"
    if arm in {"hs", "hs_stab"}:
        return "Horseshoe"
    if arm in {"cbgs_prune", "dcpf_prune"}:
        return "Redundancy pruning"
    if arm == "cbgs_adapt":
        return "Adaptive weighting"
    if arm in {"ens_rank", "ens_rank3", "kswitch", "ens2", "ens3"}:
        return "Rank ensemble"
    return "Other"


def finite(value):
    return value is not None and math.isfinite(value)


def as_float(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result


def mean(values):
    values = [x for x in values if finite(x)]
    return sum(values) / len(values) if values else float("nan")


def sample_sd(values):
    values = [x for x in values if finite(x)]
    return statistics.stdev(values) if len(values) >= 2 else float("nan")


def quantile(values, probability):
    values = sorted(x for x in values if finite(x))
    if not values:
        return float("nan")
    if len(values) == 1:
        return values[0]
    position = (len(values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return values[lower]
    weight = position - lower
    return values[lower] * (1 - weight) + values[upper] * weight


def average_ranks(values, reverse=False):
    order = sorted(range(len(values)), key=lambda i: values[i], reverse=reverse)
    ranks = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position + 1
        while end < len(order) and values[order[end]] == values[order[position]]:
            end += 1
        rank_value = (position + 1 + end) / 2
        for index in order[position:end]:
            ranks[index] = rank_value
        position = end
    return ranks


def pearson(x_values, y_values):
    pairs = [(x, y) for x, y in zip(x_values, y_values) if finite(x) and finite(y)]
    if len(pairs) < 3:
        return float("nan")
    x_values, y_values = zip(*pairs)
    mx, my = mean(x_values), mean(y_values)
    numerator = sum((x - mx) * (y - my) for x, y in pairs)
    denominator = math.sqrt(
        sum((x - mx) ** 2 for x in x_values)
        * sum((y - my) ** 2 for y in y_values)
    )
    return numerator / denominator if denominator > 0 else float("nan")


def spearman(x_values, y_values):
    pairs = [(x, y) for x, y in zip(x_values, y_values) if finite(x) and finite(y)]
    if len(pairs) < 3:
        return float("nan")
    x_values, y_values = zip(*pairs)
    return pearson(average_ranks(x_values), average_ranks(y_values))


def require_version(directory, name):
    path = directory / name
    if not path.exists():
        raise RuntimeError(f"missing provenance file: {path}")
    version = path.read_text(encoding="utf-8").strip()
    if version != CACHE_VERSION:
        raise RuntimeError(f"incompatible provenance version: {path}")


def dataset_directories(root):
    result = []
    for branch in ("validation_benchmark", "full_recipe"):
        parent = root / branch
        if parent.exists():
            result.extend((path.name, path) for path in sorted(parent.iterdir()) if path.is_dir())
    return result


def load_evaluations(name, directory, strict, evaluation_file=None):
    require_version(directory, "run_manifest.version")
    # The benchmark prediction-first table repeats the matched Random
    # rows from eval_results.csv. Load the base table first and accept only an
    # exact duplicate; any conflicting value remains a hard error.
    if evaluation_file:
        path = directory / evaluation_file
        if not path.exists():
            raise RuntimeError(f"{name}: missing deterministic evaluation: {path}")
        version_path = directory / "eval_deterministic.version"
        if not version_path.exists() or version_path.read_text(encoding="utf-8").strip() != DETERMINISTIC_EVALUATION_VERSION:
            raise RuntimeError(f"{name}: deterministic evaluation provenance is missing or incompatible")
        paths = [path]
    else:
        paths = sorted(
            (
                path for path in directory.glob("eval*.csv")
                if path.name not in {
                    "eval_deterministic.csv",
                    "eval_deterministic_manifest.csv",
                }
            ),
            key=lambda path: (path.name != "eval_results.csv", path.name),
        )
    if strict and not evaluation_file:
        observed = {path.name for path in paths}
        missing = sorted(EXPECTED_EVALUATIONS - observed)
        if missing:
            raise RuntimeError(f"{name}: missing evaluation files: {', '.join(missing)}")
    values = {}
    source = {}
    extras = {}
    for path in paths:
        if not evaluation_file and path.name != "eval_results.csv":
            extension = path.stem.removeprefix("eval_")
            require_version(directory, f"extension_manifest_{extension}.version")
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                required = {"repeat_idx", "fold_idx", "arm", "k", "AUC"}
                if not required <= row.keys():
                    continue
                key = (
                    int(row["repeat_idx"]),
                    int(row["fold_idx"]),
                    normalise_arm(row["arm"]),
                    int(row["k"]),
                )
                if key[2] in EXCLUDED_METHODS:
                    continue
                auc = as_float(row["AUC"])
                if key in values:
                    if values[key] != auc:
                        raise RuntimeError(
                            f"{name}: conflicting evaluation cell {key}"
                        )
                    continue
                if not finite(auc) or not 0 <= auc <= 1:
                    raise RuntimeError(f"{name}: invalid AUC for {key}: {row['AUC']}")
                values[key] = auc
                source[key] = path.name
                extras[key] = row
    if not values:
        raise RuntimeError(f"{name}: no evaluation rows")
    grouped = defaultdict(list)
    for (repeat_idx, fold_idx, arm, k), auc in values.items():
        grouped[(arm, k)].append(auc)
    if strict:
        methods = {arm for arm, _ in grouped if arm != "Random"}
        if len(methods) != EXPECTED_METHODS:
            raise RuntimeError(f"{name}: expected {EXPECTED_METHODS} methods, found {len(methods)}")
        for cell, observations in grouped.items():
            if len(observations) != EXPECTED_SPLITS:
                raise RuntimeError(f"{name}: incomplete evaluation {cell}: {len(observations)}")
    return values, source, extras


def predictive_summaries(dataset, values, source):
    random = {
        (repeat_idx, fold_idx, k): auc
        for (repeat_idx, fold_idx, arm, k), auc in values.items()
        if arm == "Random"
    }
    if not random:
        raise RuntimeError(f"{dataset}: matched Random rows are missing")
    methods = sorted({arm for (_, _, arm, _) in values if arm != "Random"})
    by_k_rows = []
    primary_rows = []
    sensitivity_rows = []
    method_deltas = {}

    for arm in methods:
        primary_deltas = []
        primary_aucs = []
        repeat_deltas = defaultdict(list)
        delta_by_k = {}
        for k in PANEL_SIZES:
            observations = []
            deltas = []
            files = set()
            for (repeat_idx, fold_idx, method, panel_k), auc in values.items():
                if method != arm or panel_k != k:
                    continue
                baseline = random.get((repeat_idx, fold_idx, k))
                if baseline is None:
                    raise RuntimeError(f"{dataset}: Random missing for {(repeat_idx, fold_idx, k)}")
                observations.append(auc)
                deltas.append(auc - baseline)
                files.add(source[(repeat_idx, fold_idx, method, panel_k)])
                if k in PRIMARY_K:
                    primary_deltas.append(auc - baseline)
                    primary_aucs.append(auc)
                    repeat_deltas[repeat_idx].append(auc - baseline)
            if not observations:
                continue
            delta_by_k[k] = mean(deltas)
            by_k_rows.append({
                "dataset": dataset,
                "dataset_label": LABELS.get(dataset, dataset),
                "arm": arm,
                "family": method_family(arm),
                "k": k,
                "n": len(observations),
                "mean_auc": mean(observations),
                "median_auc": quantile(observations, 0.5),
                "sd_auc": sample_sd(observations),
                "q1_auc": quantile(observations, 0.25),
                "q3_auc": quantile(observations, 0.75),
                "mean_random_auc": mean([random[key] for key in random if key[2] == k]),
                "mean_delta": mean(deltas),
                "median_delta": quantile(deltas, 0.5),
                "sd_delta": sample_sd(deltas),
                "q1_delta": quantile(deltas, 0.25),
                "q3_delta": quantile(deltas, 0.75),
                "positive_split_fraction": mean([float(value > 0) for value in deltas]),
                "source_files": ";".join(sorted(files)),
            })
        if not primary_deltas:
            continue
        repeat_means = [mean(repeat_deltas[index]) for index in sorted(repeat_deltas)]
        method_deltas[arm] = primary_deltas
        primary_rows.append({
            "dataset": dataset,
            "dataset_label": LABELS.get(dataset, dataset),
            "arm": arm,
            "family": method_family(arm),
            "primary_mean_auc": mean(primary_aucs),
            "primary_delta": mean(primary_deltas),
            "primary_median_delta": quantile(primary_deltas, 0.5),
            "primary_sd_delta": sample_sd(primary_deltas),
            "primary_q1_delta": quantile(primary_deltas, 0.25),
            "primary_q3_delta": quantile(primary_deltas, 0.75),
            "positive_split_fraction": mean([float(value > 0) for value in primary_deltas]),
            "repeat_mean_min": min(repeat_means),
            "repeat_mean_max": max(repeat_means),
            "n_split_k_cells": len(primary_deltas),
        })
        ks = sorted(delta_by_k)
        deltas = [delta_by_k[k] for k in ks]
        best_index = max(range(len(ks)), key=lambda index: deltas[index])
        worst_index = min(range(len(ks)), key=lambda index: deltas[index])
        sensitivity_rows.append({
            "dataset": dataset,
            "arm": arm,
            "family": method_family(arm),
            "best_k": ks[best_index],
            "best_delta": deltas[best_index],
            "worst_k": ks[worst_index],
            "worst_delta": deltas[worst_index],
            "delta_range": max(deltas) - min(deltas),
            "spearman_log2_k_delta": spearman([math.log2(k) for k in ks], deltas),
        })

    rank_values = [row["primary_delta"] for row in primary_rows]
    ranks = average_ranks(rank_values, reverse=True)
    for row, rank in zip(primary_rows, ranks):
        row["dataset_rank"] = rank

    pairwise_rows = []
    for arm_a, arm_b in itertools.combinations(methods, 2):
        differences = []
        repeat_differences = defaultdict(list)
        for repeat_idx in (1, 2, 3):
            for fold_idx in (1, 2, 3, 4, 5):
                for k in PRIMARY_K:
                    a = values.get((repeat_idx, fold_idx, arm_a, k))
                    b = values.get((repeat_idx, fold_idx, arm_b, k))
                    if a is None or b is None:
                        continue
                    differences.append(a - b)
                    repeat_differences[repeat_idx].append(a - b)
        if not differences:
            continue
        repeat_means = [mean(repeat_differences[index]) for index in sorted(repeat_differences)]
        pairwise_rows.append({
            "dataset": dataset,
            "arm_a": arm_a,
            "arm_b": arm_b,
            "mean_auc_difference_a_minus_b": mean(differences),
            "median_difference": quantile(differences, 0.5),
            "sd_difference": sample_sd(differences),
            "q1_difference": quantile(differences, 0.25),
            "q3_difference": quantile(differences, 0.75),
            "a_win_fraction": mean([float(value > 0) for value in differences]),
            "tie_fraction": mean([float(value == 0) for value in differences]),
            "repeat_mean_min": min(repeat_means),
            "repeat_mean_max": max(repeat_means),
            "n_split_k_cells": len(differences),
        })
    return by_k_rows, primary_rows, sensitivity_rows, pairwise_rows


def read_ranking(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return [row["gene"] for row in rows if row.get("gene")]


def stability_summaries(dataset, directory, methods, strict):
    rows = []
    ungrouped = "GS_full_ungrouped" if directory.parent.name == "validation_benchmark" else "full_ungrouped"
    for arm in methods:
        rankings = []
        for repeat_idx in (1, 2, 3):
            for fold_idx in (1, 2, 3, 4, 5):
                raw_arm = arm
                if directory.parent.name == "full_recipe":
                    raw_arm = {"GS_full_grouped": "full_grouped", "GS_full_ungrouped": "full_ungrouped"}.get(arm, arm)
                path = directory / f"ranking_r{repeat_idx}_f{fold_idx}_{raw_arm}.csv"
                if not path.exists():
                    if strict:
                        raise RuntimeError(f"{dataset}: missing ranking {path.name}")
                    rankings = []
                    break
                ranking = read_ranking(path)
                if not ranking and raw_arm in {"GS_full_grouped", "full_grouped"}:
                    fallback = directory / f"ranking_r{repeat_idx}_f{fold_idx}_{ungrouped}.csv"
                    ranking = read_ranking(fallback)
                if len(ranking) != len(set(ranking)):
                    raise RuntimeError(f"{dataset}: duplicate genes in {path.name}")
                rankings.append(ranking)
            if not rankings:
                break
        if len(rankings) != EXPECTED_SPLITS:
            continue
        universe = set().union(*(set(ranking) for ranking in rankings))
        d = len(universe)
        for k in PANEL_SIZES:
            panels = [set(ranking[:k]) for ranking in rankings]
            if any(len(panel) == 0 for panel in panels):
                continue
            intersections = []
            jaccards = []
            for panel_a, panel_b in itertools.combinations(panels, 2):
                intersection = len(panel_a & panel_b)
                union = len(panel_a | panel_b)
                intersections.append(intersection)
                jaccards.append(intersection / union if union else float("nan"))
            counts = Counter(gene for panel in panels for gene in panel)
            probabilities = [counts.get(gene, 0) / len(panels) for gene in universe]
            sample_variances = [
                len(panels) / (len(panels) - 1) * probability * (1 - probability)
                for probability in probabilities
            ]
            q = mean([len(panel) for panel in panels])
            denominator = (q / d) * (1 - q / d) if d > 0 else 0
            nogueira = 1 - mean(sample_variances) / denominator if denominator > 0 else float("nan")
            rows.append({
                "dataset": dataset,
                "arm": arm,
                "family": method_family(arm),
                "k": k,
                "n_rankings": len(rankings),
                "ranking_universe": d,
                "mean_pairwise_jaccard": mean(jaccards),
                "median_pairwise_jaccard": quantile(jaccards, 0.5),
                "q1_pairwise_jaccard": quantile(jaccards, 0.25),
                "q3_pairwise_jaccard": quantile(jaccards, 0.75),
                "mean_pairwise_overlap": mean(intersections),
                "nogueira_stability": nogueira,
                "genes_in_all_splits": sum(count == len(panels) for count in counts.values()),
                "genes_in_at_least_half_splits": sum(count >= math.ceil(len(panels) / 2) for count in counts.values()),
            })
    return rows


def runtime_summaries(dataset, directory):
    pattern = re.compile(r"ranking_r(\d+)_f(\d+)_(.+)_meta[.]csv$")
    grouped = defaultdict(list)
    for path in sorted(directory.glob("ranking_r*_f*_*_meta.csv")):
        match = pattern.match(path.name)
        if not match:
            continue
        arm = normalise_arm(match.group(3))
        with path.open(newline="", encoding="utf-8") as handle:
            metadata = {row["key"]: row["value"] for row in csv.DictReader(handle)}
        seconds = as_float(metadata.get("fit_seconds"))
        if finite(seconds):
            grouped[arm].append(seconds)
    return [{
        "dataset": dataset,
        "arm": arm,
        "family": method_family(arm),
        "n_fits": len(seconds),
        "mean_fit_seconds": mean(seconds),
        "median_fit_seconds": quantile(seconds, 0.5),
        "q1_fit_seconds": quantile(seconds, 0.25),
        "q3_fit_seconds": quantile(seconds, 0.75),
        "total_fit_seconds": sum(seconds),
    } for arm, seconds in sorted(grouped.items())]


def gate_summaries(dataset, values, extras):
    observed = {}
    for key, row in extras.items():
        arm = key[2]
        if arm != "GS_full_grouped":
            continue
        value = row.get("gate_fallback", "")
        if value != "":
            observed[(key[0], key[1], arm)] = value.lower() == "true"
    grouped = defaultdict(list)
    for (_, _, arm), flag in observed.items():
        grouped[arm].append(flag)
    return [{
        "dataset": dataset,
        "arm": arm,
        "n_splits": len(flags),
        "fallback_fraction": mean([float(flag) for flag in flags]),
    } for arm, flags in sorted(grouped.items())]


def adaptive_weight_summaries(dataset, extras):
    observed = {}
    for (repeat_idx, fold_idx, arm, _), row in extras.items():
        if arm != "cbgs_adapt" or not row.get("w_chosen"):
            continue
        observed[(repeat_idx, fold_idx)] = as_float(row["w_chosen"])
    counts = Counter(observed.values())
    return [{
        "dataset": dataset,
        "w_chosen": weight,
        "split_count": count,
        "split_fraction": count / len(observed),
    } for weight, count in sorted(counts.items())] if observed else []


def biology_summaries(dataset, directory, prediction_by_k, strict):
    path = directory / "biology_string.csv"
    if not path.exists():
        if strict:
            raise RuntimeError(f"{dataset}: biology_string.csv is missing")
        return [], []
    require_version(directory, "extension_manifest_biology.version")
    grouped = defaultdict(list)
    keys = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            arm = normalise_arm(row["arm"])
            if arm in EXCLUDED_METHODS:
                continue
            key = (int(row["repeat_idx"]), int(row["fold_idx"]), arm)
            if key in keys:
                raise RuntimeError(f"{dataset}: duplicate biology cell {key}")
            keys.add(key)
            grouped[arm].append(row)
    if strict and len(grouped) != EXPECTED_METHODS:
        raise RuntimeError(f"{dataset}: expected biology for {EXPECTED_METHODS} methods, found {len(grouped)}")
    summary = []
    for arm, observations in sorted(grouped.items()):
        if strict and len(observations) != EXPECTED_SPLITS:
            raise RuntimeError(f"{dataset}: incomplete biology for {arm}: {len(observations)}")
        ratios = [as_float(row.get("ratio")) for row in observations]
        observed_edges = [as_float(row.get("obs_edges")) for row in observations]
        expected_edges = [as_float(row.get("exp_edges")) for row in observations]
        mapped = [as_float(row.get("n_mapped")) for row in observations]
        panel_sizes = [as_float(row.get("n_panel", 50)) for row in observations]
        prediction = prediction_by_k.get((dataset, arm, 50), {})
        summary.append({
            "dataset": dataset,
            "arm": arm,
            "family": method_family(arm),
            "n_splits": len(observations),
            "mean_string_ratio": mean(ratios),
            "median_string_ratio": quantile(ratios, 0.5),
            "sd_string_ratio": sample_sd(ratios),
            "q1_string_ratio": quantile(ratios, 0.25),
            "q3_string_ratio": quantile(ratios, 0.75),
            "mean_observed_edges": mean(observed_edges),
            "mean_expected_edges": mean(expected_edges),
            "mean_mapped_genes": mean(mapped),
            "mean_mapping_fraction": mean([
                mapped_value / panel_size
                for mapped_value, panel_size in zip(mapped, panel_sizes)
                if finite(mapped_value) and finite(panel_size) and panel_size > 0
            ]),
            "k50_mean_auc": prediction.get("mean_auc", float("nan")),
            "k50_delta": prediction.get("mean_delta", float("nan")),
        })
    biology_ranks = average_ranks([row["mean_string_ratio"] for row in summary], reverse=True)
    prediction_ranks = average_ranks([row["k50_delta"] for row in summary], reverse=True)
    for row, biology_rank, prediction_rank in zip(summary, biology_ranks, prediction_ranks):
        row["biology_rank"] = biology_rank
        row["prediction_rank_k50"] = prediction_rank
        row["rank_difference_biology_minus_prediction"] = biology_rank - prediction_rank
    association = [{
        "dataset": dataset,
        "n_methods": len(summary),
        "spearman_string_ratio_vs_k50_delta": spearman(
            [row["mean_string_ratio"] for row in summary],
            [row["k50_delta"] for row in summary],
        ),
        "pearson_string_ratio_vs_k50_delta": pearson(
            [row["mean_string_ratio"] for row in summary],
            [row["k50_delta"] for row in summary],
        ),
    }]
    return summary, association


def disease_biology_summaries(dataset, directory, prediction_by_k, strict):
    path = directory / "biology_disease.csv"
    if not path.exists():
        if strict:
            raise RuntimeError(f"{dataset}: biology_disease.csv is missing")
        return [], []
    require_version(directory, "extension_manifest_biology_disease.version")
    grouped = defaultdict(list)
    keys = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            arm = normalise_arm(row["arm"])
            if arm in EXCLUDED_METHODS:
                continue
            key = (int(row["repeat_idx"]), int(row["fold_idx"]), arm)
            if key in keys:
                raise RuntimeError(f"{dataset}: duplicate disease-biology cell {key}")
            keys.add(key)
            grouped[arm].append(row)
    if strict and len(grouped) != EXPECTED_METHODS:
        raise RuntimeError(
            f"{dataset}: expected disease biology for {EXPECTED_METHODS} methods, "
            f"found {len(grouped)}"
        )
    summary = []
    for arm, observations in sorted(grouped.items()):
        if strict and len(observations) != EXPECTED_SPLITS:
            raise RuntimeError(
                f"{dataset}: incomplete disease biology for {arm}: {len(observations)}"
            )
        ratios = [as_float(row.get("enrichment_ratio")) for row in observations]
        z_scores = [as_float(row.get("enrichment_z")) for row in observations]
        overlaps = [as_float(row.get("n_target_overlap")) for row in observations]
        recalls = [as_float(row.get("target_recall")) for row in observations]
        empirical = [as_float(row.get("empirical_p")) for row in observations]
        prediction = prediction_by_k.get((dataset, arm, 50), {})
        summary.append({
            "dataset": dataset,
            "arm": arm,
            "family": method_family(arm),
            "disease_label": observations[0].get("disease_label", ""),
            "ontology_id": observations[0].get("ontology_id", ""),
            "n_splits": len(observations),
            "mean_disease_enrichment_ratio": mean(ratios),
            "median_disease_enrichment_ratio": quantile(ratios, 0.5),
            "q1_disease_enrichment_ratio": quantile(ratios, 0.25),
            "q3_disease_enrichment_ratio": quantile(ratios, 0.75),
            "mean_disease_enrichment_z": mean(z_scores),
            "median_disease_enrichment_z": quantile(z_scores, 0.5),
            "mean_target_overlap": mean(overlaps),
            "mean_target_recall": mean(recalls),
            "median_empirical_p": quantile(empirical, 0.5),
            "k50_mean_auc": prediction.get("mean_auc", float("nan")),
            "k50_delta": prediction.get("mean_delta", float("nan")),
        })
    disease_ranks = average_ranks(
        [row["mean_disease_enrichment_z"] for row in summary], reverse=True
    )
    prediction_ranks = average_ranks(
        [row["k50_delta"] for row in summary], reverse=True
    )
    for row, disease_rank, prediction_rank in zip(summary, disease_ranks, prediction_ranks):
        row["disease_biology_rank"] = disease_rank
        row["prediction_rank_k50"] = prediction_rank
        row["disease_rank_minus_prediction_rank"] = disease_rank - prediction_rank
    association = [{
        "dataset": dataset,
        "n_disease_methods": len(summary),
        "spearman_disease_z_vs_k50_delta": spearman(
            [row["mean_disease_enrichment_z"] for row in summary],
            [row["k50_delta"] for row in summary],
        ),
        "pearson_disease_z_vs_k50_delta": pearson(
            [row["mean_disease_enrichment_z"] for row in summary],
            [row["k50_delta"] for row in summary],
        ),
    }]
    return summary, association


def overall_method_summaries(primary_rows):
    datasets = sorted({row["dataset"] for row in primary_rows})
    grouped = defaultdict(list)
    for row in primary_rows:
        grouped[row["arm"]].append(row)
    result = []
    for arm, rows in sorted(grouped.items()):
        ranks = [row["dataset_rank"] for row in rows]
        deltas = [row["primary_delta"] for row in rows]
        result.append({
            "arm": arm,
            "family": method_family(arm),
            "datasets_present": len(rows),
            "datasets_total": len(datasets),
            "mean_dataset_rank": mean(ranks),
            "median_dataset_rank": quantile(ranks, 0.5),
            "first_place_count": sum(rank == 1 for rank in ranks),
            "mean_primary_delta_across_datasets": mean(deltas),
            "median_primary_delta_across_datasets": quantile(deltas, 0.5),
            "minimum_primary_delta": min(deltas),
            "maximum_primary_delta": max(deltas),
        })
    return sorted(result, key=lambda row: (row["mean_dataset_rank"], -row["mean_primary_delta_across_datasets"]))


def write_csv(path, rows):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def fmt(value, digits=3, signed=False):
    if not finite(value):
        return "NA"
    return f"{value:+.{digits}f}" if signed else f"{value:.{digits}f}"


def render_markdown(path, primary_rows, overall_rows, biology_rows, associations, stability_rows, strict, deterministic):
    lines = [
        "# GeneSelectR redesign comprehensive comparison",
        "",
        "## Summary",
        "",
        "Primary performance is the arithmetic mean of AUC minus matched Random at k = 10, 20, and 50. "
        "The results are descriptive because repeated cross-validation folds share observations.",
        "",
        f"Report status: {'complete deterministic benchmark' if deterministic and strict else 'complete corrected benchmark' if strict else 'interim corrected benchmark'}.",
        "",
        "The ensemble used explicit stratified glmnet folds and explicit seeds for glmnet, XGBoost, and ranger. "
        "The dc_pf_ctrl duplicate was retained as a file-level reproducibility check and excluded from method rankings.",
        "",
        "## Cross-dataset comparison",
        "",
        "| method | datasets | mean rank | first-place datasets | mean primary delta | minimum delta |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in overall_rows:
        lines.append(
            f"| {row['arm']} | {row['datasets_present']} | {fmt(row['mean_dataset_rank'], 2)} | "
            f"{row['first_place_count']} | {fmt(row['mean_primary_delta_across_datasets'], signed=True)} | "
            f"{fmt(row['minimum_primary_delta'], signed=True)} |"
        )
    for dataset in sorted({row["dataset"] for row in primary_rows}):
        rows = sorted(
            [row for row in primary_rows if row["dataset"] == dataset],
            key=lambda row: row["dataset_rank"],
        )
        lines.extend([
            "",
            f"## {LABELS.get(dataset, dataset)}",
            "",
            "| rank | method | primary delta | median | IQR | positive cells | repeat range |",
            "|---:|---|---:|---:|---:|---:|---:|",
        ])
        for row in rows:
            lines.append(
                f"| {fmt(row['dataset_rank'], 1)} | {row['arm']} | {fmt(row['primary_delta'], signed=True)} | "
                f"{fmt(row['primary_median_delta'], signed=True)} | "
                f"{fmt(row['primary_q1_delta'], signed=True)} to {fmt(row['primary_q3_delta'], signed=True)} | "
                f"{100 * row['positive_split_fraction']:.1f}% | "
                f"{fmt(row['repeat_mean_min'], signed=True)} to {fmt(row['repeat_mean_max'], signed=True)} |"
            )
        bio = sorted(
            [row for row in biology_rows if row["dataset"] == dataset],
            key=lambda row: row["biology_rank"],
        )
        if bio:
            association = next(row for row in associations if row["dataset"] == dataset)
            lines.extend([
                "",
                "### Biology",
                "",
                "STRING connectivity was calculated for each exact top-50 panel against 1,000 random panels from the same split-specific variance pool. "
                "Disease-specific enrichment used frozen Open Targets associations and the same candidate-pool null. "
                f"The method-level Spearman correlations with predictive delta were {fmt(association['spearman_string_ratio_vs_k50_delta'])} for STRING and "
                f"{fmt(association.get('spearman_disease_z_vs_k50_delta', float('nan')))} for disease-association enrichment.",
                "",
                "| STRING rank | method | STRING ratio | disease z | target overlap | mapped genes | k=50 predictive delta |",
                "|---:|---|---:|---:|---:|---:|---:|",
            ])
            for row in bio:
                lines.append(
                    f"| {fmt(row['biology_rank'], 1)} | {row['arm']} | {fmt(row['mean_string_ratio'], 2)} | "
                    f"{fmt(row.get('mean_disease_enrichment_z', float('nan')), 2)} | "
                    f"{fmt(row.get('mean_target_overlap', float('nan')), 1)} | "
                    f"{fmt(row['mean_mapped_genes'], 1)} | {fmt(row['k50_delta'], signed=True)} |"
                )
    lines.extend([
        "",
        "## Stability",
        "",
        "Stability tables are stored in `stability_comparisons.csv`. Nogueira stability and pairwise Jaccard overlap use the 15 top-k panels. "
        "Training sets overlap across folds, so these estimates are descriptive and are expected to be inflated.",
        "",
        "## Interpretation limits",
        "",
        "- SOS-ALL was used for method development. The other six datasets were used for validation after development.",
        "- STRING connectivity measures network coherence. Open Targets enrichment supplies a separate disease-specific axis. Neither biology measure is combined with predictive AUC.",
        "- Disease-specific mappings were specified after predictive results were available and are exploratory.",
        "- Pairwise differences are reported without p-values because the repeated cross-validation cells are dependent.",
        "- The psoriasis positive control is subject to an AUC ceiling.",
        "- SOS-ALL is retained as a negative-control dataset.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


HTML_TEMPLATE = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GeneSelectR redesign comparison</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.5.1"></script>
<style>
:root{--bg:#f4f6f8;--card:#fff;--ink:#17212b;--muted:#5f6b76;--line:#dce2e7;--blue:#2f6f9f;--orange:#c66a2b;--green:#3e8061;--red:#a84b4b}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
.wrap{max-width:1500px;margin:auto;padding:18px}.header{background:#183246;color:#fff;padding:20px 24px;border-radius:9px;display:flex;gap:20px;justify-content:space-between;align-items:flex-end;flex-wrap:wrap}
h1{font-size:22px;margin:0 0 5px}h2{font-size:17px;margin:0 0 14px}.sub{color:#d8e2ea}.filters{display:flex;gap:12px;flex-wrap:wrap}.filter label{display:block;font-size:11px;color:#c9d5de;margin-bottom:3px}.filter select{min-width:190px;padding:7px 9px;border:1px solid #6d8291;border-radius:5px;background:#fff;color:#17212b}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px;margin:14px 0}.card,.panel{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:16px}.klabel{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}.kvalue{font-size:27px;font-weight:700;margin-top:4px}
.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin-bottom:14px}.chart{height:380px}.wide{grid-column:1/-1}.note{color:var(--muted);font-size:12px;margin-top:7px}.tables{display:grid;gap:14px}.scroll{overflow:auto;max-height:600px}table{width:100%;border-collapse:collapse;font-size:12px}th{position:sticky;top:0;background:#eef2f5;text-align:left;padding:8px;border-bottom:2px solid var(--line);white-space:nowrap;cursor:pointer}td{padding:7px 8px;border-bottom:1px solid #edf0f2;white-space:nowrap}tr:hover td{background:#f7f9fa}.positive{color:var(--green)}.negative{color:var(--red)}
.status{display:inline-block;padding:3px 7px;border-radius:10px;background:#dceaf3;color:#214c6a;font-size:11px;font-weight:600}.methods{font-size:12px;color:var(--muted);margin-top:8px}
@media(max-width:850px){.grid{grid-template-columns:1fr}.wide{grid-column:auto}.chart{height:320px}}@media print{body{background:#fff}.filters{display:none}.panel,.card{break-inside:avoid}}
</style></head><body><div class="wrap">
<div class="header"><div><h1>GeneSelectR redesign comparison</h1><div class="sub">Prediction, panel-size sensitivity, stability, runtime, STRING connectivity, and disease associations</div></div><div class="filters"><div class="filter"><label>Dataset</label><select id="dataset"></select></div><div class="filter"><label>Method family</label><select id="family"><option value="all">All families</option></select></div></div></div>
<div class="kpis"><div class="card"><div class="klabel">Datasets</div><div class="kvalue" id="nDataset"></div></div><div class="card"><div class="klabel">Methods in selected dataset</div><div class="kvalue" id="nMethod"></div></div><div class="card"><div class="klabel">Evaluation cells</div><div class="kvalue" id="nCells"></div></div><div class="card"><div class="klabel">Biology panels</div><div class="kvalue" id="nBio"></div></div></div>
<div class="grid"><div class="panel chart"><h2>Primary AUC minus matched Random</h2><canvas id="primaryChart"></canvas><div class="note">Mean across k=10, 20, and 50. Bars show all methods in the selected family.</div></div><div class="panel chart"><h2>Panel-size sensitivity</h2><canvas id="sizeChart"></canvas><div class="note">The eight highest primary-ranked methods are displayed.</div></div><div class="panel chart"><h2>STRING connectivity and prediction at k=50</h2><canvas id="biologyChart"></canvas><div class="note">STRING ratio uses 1,000 candidate-pool-matched random panels per split.</div></div><div class="panel chart"><h2>Disease association and prediction at k=50</h2><canvas id="diseaseChart"></canvas><div class="note">Disease enrichment uses frozen Open Targets associations and the same candidate-pool null.</div></div></div>
<div class="tables"><div class="panel"><h2>Predictive comparison</h2><div class="scroll"><table id="primaryTable"></table></div></div><div class="panel"><h2>Biology comparison</h2><div class="scroll"><table id="biologyTable"></table></div></div><div class="panel"><h2>Methodology and limits</h2><span class="status">__STATUS__</span><p>Repeated cross-validation cells share observations. Means, medians, interquartile ranges, repeat ranges, and paired win fractions are descriptive. Hypothesis-test p-values are omitted.</p><p>SOS-ALL was used for method development. The other six datasets were used for validation after development. STRING measures network coherence. Open Targets provides a separate disease-specific association axis. The disease mappings were specified after predictive results were available and are exploratory.</p><div class="methods">Detailed machine-readable files include every gene-set size, every pairwise method comparison, stability, runtime, gate fallback, adaptive weights, and biology-prediction associations.</div></div></div>
</div><script>
const DATA=__DATA__;
const COLORS=['#2f6f9f','#c66a2b','#3e8061','#8f5a9d','#a84b4b','#607d8b','#9a7b34','#477d86'];
const dsSel=document.getElementById('dataset'),famSel=document.getElementById('family');
const datasets=[...new Set(DATA.primary.map(d=>d.dataset))].sort();datasets.forEach(d=>{const o=document.createElement('option');o.value=d;o.textContent=(DATA.labels[d]||d);dsSel.appendChild(o)});
const families=[...new Set(DATA.primary.map(d=>d.family))].sort();families.forEach(f=>{const o=document.createElement('option');o.value=f;o.textContent=f;famSel.appendChild(o)});
let charts={};function f3(x,sign=false){if(x===null||!Number.isFinite(x))return 'NA';return (sign&&x>=0?'+':'')+x.toFixed(3)}
function destroy(name){if(charts[name])charts[name].destroy()}
function filtered(rows){return rows.filter(r=>r.dataset===dsSel.value&&(famSel.value==='all'||r.family===famSel.value))}
function table(id,cols,rows){const el=document.getElementById(id);let sort=cols[0][0],dir=1;function draw(){const data=[...rows].sort((a,b)=>{const x=a[sort],y=b[sort];return dir*((x??'')<(y??'')?-1:(x??'')>(y??'')?1:0)});el.innerHTML='<thead><tr>'+cols.map(c=>`<th data-k="${c[0]}">${c[1]}</th>`).join('')+'</tr></thead><tbody>'+data.map(r=>'<tr>'+cols.map(c=>`<td>${c[2]?c[2](r[c[0]]):r[c[0]]}</td>`).join('')+'</tr>').join('')+'</tbody>';el.querySelectorAll('th').forEach(h=>h.onclick=()=>{if(sort===h.dataset.k)dir*=-1;else{sort=h.dataset.k;dir=1}draw()})}draw()}
function render(){const ds=dsSel.value;const p=filtered(DATA.primary).sort((a,b)=>b.primary_delta-a.primary_delta);const byk=filtered(DATA.by_k);const bio=filtered(DATA.biology);document.getElementById('nDataset').textContent=datasets.length;document.getElementById('nMethod').textContent=p.length;document.getElementById('nCells').textContent=byk.reduce((s,r)=>s+r.n,0).toLocaleString();document.getElementById('nBio').textContent=bio.length;
destroy('p');charts.p=new Chart(document.getElementById('primaryChart'),{type:'bar',data:{labels:p.map(r=>r.arm),datasets:[{data:p.map(r=>r.primary_delta),backgroundColor:p.map(r=>r.primary_delta>=0?'#3e8061cc':'#a84b4bcc')}]},options:{indexAxis:'y',responsive:true,maintainAspectRatio:false,animation:false,plugins:{legend:{display:false}},scales:{x:{title:{display:true,text:'AUC minus Random'}}}}});
const top=p.slice(0,8).map(r=>r.arm);destroy('s');charts.s=new Chart(document.getElementById('sizeChart'),{type:'line',data:{labels:[10,20,50,100,200,500],datasets:top.map((m,i)=>({label:m,data:[10,20,50,100,200,500].map(k=>{const r=byk.find(x=>x.arm===m&&x.k===k);return r?r.mean_delta:null}),borderColor:COLORS[i%COLORS.length],backgroundColor:COLORS[i%COLORS.length],tension:.15,pointRadius:3}))},options:{responsive:true,maintainAspectRatio:false,animation:false,scales:{x:{type:'logarithmic',title:{display:true,text:'Panel size k'}},y:{title:{display:true,text:'AUC minus Random'}}}}});
destroy('b');charts.b=new Chart(document.getElementById('biologyChart'),{type:'scatter',data:{datasets:bio.map((r,i)=>({label:r.arm,data:[{x:r.mean_string_ratio,y:r.k50_delta}],backgroundColor:COLORS[i%COLORS.length],pointRadius:6}))},options:{responsive:true,maintainAspectRatio:false,animation:false,plugins:{legend:{position:'right'}},scales:{x:{title:{display:true,text:'Mean STRING enrichment ratio'}},y:{title:{display:true,text:'AUC minus Random at k=50'}}}}});
destroy('d');charts.d=new Chart(document.getElementById('diseaseChart'),{type:'scatter',data:{datasets:bio.filter(r=>Number.isFinite(r.mean_disease_enrichment_z)).map((r,i)=>({label:r.arm,data:[{x:r.mean_disease_enrichment_z,y:r.k50_delta}],backgroundColor:COLORS[i%COLORS.length],pointRadius:6}))},options:{responsive:true,maintainAspectRatio:false,animation:false,plugins:{legend:{position:'right'}},scales:{x:{title:{display:true,text:'Mean disease-association enrichment z score'}},y:{title:{display:true,text:'AUC minus Random at k=50'}}}}});
table('primaryTable',[['dataset_rank','Rank',x=>Number(x).toFixed(1)],['arm','Method'],['family','Family'],['primary_delta','Primary delta',x=>f3(x,true)],['primary_median_delta','Median',x=>f3(x,true)],['primary_q1_delta','Q1',x=>f3(x,true)],['primary_q3_delta','Q3',x=>f3(x,true)],['positive_split_fraction','Positive cells',x=>(100*x).toFixed(1)+'%'],['repeat_mean_min','Repeat minimum',x=>f3(x,true)],['repeat_mean_max','Repeat maximum',x=>f3(x,true)]],p);
table('biologyTable',[['biology_rank','STRING rank',x=>Number(x).toFixed(1)],['disease_biology_rank','Disease rank',x=>x==null?'NA':Number(x).toFixed(1)],['arm','Method'],['mean_string_ratio','STRING ratio',x=>Number(x).toFixed(2)],['mean_disease_enrichment_z','Disease z',x=>x==null?'NA':Number(x).toFixed(2)],['mean_target_overlap','Target overlap',x=>x==null?'NA':Number(x).toFixed(1)],['mean_mapped_genes','STRING mapped',x=>Number(x).toFixed(1)],['k50_delta','k=50 delta',x=>f3(x,true)]],bio)}
dsSel.onchange=render;famSel.onchange=render;if(datasets.length){dsSel.value=datasets[0];render()}
</script></body></html>'''


def json_safe(rows):
    cleaned = []
    for row in rows:
        cleaned.append({key: (value if not isinstance(value, float) or finite(value) else None) for key, value in row.items()})
    return cleaned


def main():
    args = parse_args()
    root = Path(args.root)
    strict = not args.allow_partial
    output_dir = Path(args.output_dir) if args.output_dir else root / "full_benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_by_k = []
    all_primary = []
    all_sensitivity = []
    all_pairwise = []
    all_stability = []
    all_runtime = []
    all_gate = []
    all_weights = []
    dataset_records = []

    for dataset, directory in dataset_directories(root):
        if args.allow_partial and not any(directory.glob("eval*.csv")):
            continue
        values, source, extras = load_evaluations(
            dataset, directory, strict, args.evaluation_file
        )
        by_k, primary, sensitivity, pairwise = predictive_summaries(dataset, values, source)
        methods = sorted({row["arm"] for row in primary})
        all_by_k.extend(by_k)
        all_primary.extend(primary)
        all_sensitivity.extend(sensitivity)
        all_pairwise.extend(pairwise)
        all_stability.extend(stability_summaries(dataset, directory, methods, strict))
        all_runtime.extend(runtime_summaries(dataset, directory))
        all_gate.extend(gate_summaries(dataset, values, extras))
        all_weights.extend(adaptive_weight_summaries(dataset, extras))
        dataset_records.append((dataset, directory))

    if not dataset_records:
        raise RuntimeError(f"no result datasets under {root}")

    prediction_by_k = {
        (row["dataset"], row["arm"], row["k"]): row for row in all_by_k
    }
    all_string_biology = []
    all_disease_biology = []
    association_by_dataset = {}
    for dataset, directory in dataset_records:
        string_biology, string_association = biology_summaries(
            dataset, directory, prediction_by_k, strict
        )
        disease_biology, disease_association = disease_biology_summaries(
            dataset, directory, prediction_by_k, strict
        )
        all_string_biology.extend(string_biology)
        all_disease_biology.extend(disease_biology)
        association_by_dataset.setdefault(dataset, {"dataset": dataset})
        if string_association:
            association_by_dataset[dataset].update(string_association[0])
        if disease_association:
            association_by_dataset[dataset].update(disease_association[0])

    biology_by_method = {
        (row["dataset"], row["arm"]): dict(row) for row in all_string_biology
    }
    for row in all_disease_biology:
        key = (row["dataset"], row["arm"])
        if key not in biology_by_method:
            biology_by_method[key] = {
                "dataset": row["dataset"],
                "arm": row["arm"],
                "family": row["family"],
                "biology_rank": float("nan"),
                "mean_string_ratio": float("nan"),
                "q1_string_ratio": float("nan"),
                "q3_string_ratio": float("nan"),
                "mean_mapped_genes": float("nan"),
                "k50_delta": row["k50_delta"],
            }
        biology_by_method[key].update(row)
    all_biology = sorted(
        biology_by_method.values(), key=lambda row: (row["dataset"], row["arm"])
    )
    all_associations = [association_by_dataset[key] for key in sorted(association_by_dataset)]

    overall = overall_method_summaries(all_primary)
    outputs = {
        "predictive_by_panel_size.csv": all_by_k,
        "primary_comparisons.csv": all_primary,
        "all_pairwise_primary_comparisons.csv": all_pairwise,
        "panel_size_sensitivity.csv": all_sensitivity,
        "cross_dataset_method_summary.csv": overall,
        "stability_comparisons.csv": all_stability,
        "runtime_comparisons.csv": all_runtime,
        "gate_fallback_summary.csv": all_gate,
        "adaptive_weight_summary.csv": all_weights,
        "biology_string_summary.csv": all_string_biology,
        "biology_disease_summary.csv": all_disease_biology,
        "biology_prediction_association.csv": all_associations,
    }
    for filename, rows in outputs.items():
        write_csv(output_dir / filename, rows)

    render_markdown(
        output_dir / "COMPREHENSIVE_REPORT.md",
        all_primary,
        overall,
        all_biology,
        all_associations,
        all_stability,
        strict,
        bool(args.evaluation_file),
    )
    payload = {
        "labels": LABELS,
        "primary": json_safe(all_primary),
        "by_k": json_safe(all_by_k),
        "biology": json_safe(all_biology),
    }
    html = HTML_TEMPLATE.replace("__DATA__", json.dumps(payload, separators=(",", ":")))
    status = "Complete deterministic benchmark" if args.evaluation_file and strict else "Complete corrected benchmark" if strict else "Interim corrected benchmark"
    html = html.replace("__STATUS__", status)
    (output_dir / "COMPREHENSIVE_REPORT.html").write_text(html, encoding="utf-8")
    print(f"wrote {output_dir / 'COMPREHENSIVE_REPORT.md'}")
    print(f"wrote {output_dir / 'COMPREHENSIVE_REPORT.html'}")
    for filename, rows in outputs.items():
        print(f"wrote {output_dir / filename} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
