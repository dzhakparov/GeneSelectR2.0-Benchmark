#!/usr/bin/env python3
"""
Fetch and standardise GEO datasets for the GeneSelectR benchmark.

For each accession this produces three files in <outdir>/<accession>/:

    expression.tsv   genes x samples, gene SYMBOLS as row names
    metadata.tsv     samples x attributes, parsed from !Sample_characteristics
    manifest.json    platform, scale, dimensions, and the checks that were run

The manifest's `scale` field is the one the R side branches on:

    "counts"            -> integer counts; needs edgeR TMM + log-CPM
    "log_normalized"    -> already log-scale; use as-is + variance filter
    "linear_normalized" -> normalized but NOT logged; needs log2(x + 1)
    "linear_background_corrected" -> linear array intensity with negative
                                      background estimates; floor at zero,
                                      then log2(x + 1)

That field is DETECTED from the matrix, not taken on trust, because a file
named "counts" containing TPM is a silent disaster: TMM normalisation of TPM
is wrong rather than merely suboptimal, and nothing downstream would complain.

WHY PROBE MAPPING MATTERS
-------------------------
Microarray series matrices are indexed by platform probe IDs ("1007_s_at"), not
gene symbols. The benchmark's biology pillar looks genes up in GO, STRING and
Open Targets by symbol, so an unmapped matrix silently scores zero biology for
every gene. This script fetches the GPL annotation and collapses probes to
symbols, keeping the highest-mean probe per symbol (the standard choice: the
mean across probes is dragged down by non-responsive probes for the same gene).

USAGE
-----
    python fetch_geo_data.py --list
    python fetch_geo_data.py GSE65682
    python fetch_geo_data.py --all --outdir data/
    python fetch_geo_data.py GSE65682 --inspect

NOT COVERED
-----------
IMvigor210 (easierData / GitHub RDS, handled in R) and the Liu/Gide melanoma
cohorts (TPM-only as supplementary material; raw counts are in dbGaP under
controlled access). Those are not GEO series matrices and are out of scope here.
"""

import argparse
import gzip
import io
import json
import re
import shutil
import sys
import urllib.error
import urllib.request
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("This script needs pandas:  pip install pandas")


GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo"
USER_AGENT = "GeneSelectR-benchmark-fetcher/1.0"


# =============================================================================
#  Dataset registry
# =============================================================================
#
#  `expected_scale` records what the literature says the series ships as, and is
#  used ONLY to cross-check the detected scale. A mismatch is reported loudly
#  rather than silently accepted.
#
#  `outcome_hint` is a regex matched against the parsed characteristic COLUMN
#  NAMES to suggest which one holds the class label. It is a hint, not a
#  decision: the script prints every characteristic it found so the label can be
#  chosen by eye. Guessing the outcome column silently is how you end up
#  benchmarking against the wrong variable.

DATASETS = {
    "GSE65682": dict(
        label="Sepsis, MARS cohort (whole blood)",
        platform="GPL13667",              # Affymetrix Human Genome U219
        expected_scale="log_normalized",
        outcome_hint=r"mortality|survival|day|outcome|death",
        notes=("802 samples: 42 healthy controls + 760 sepsis. The useful "
               "contrast is 28-day survivor vs non-survivor (~114 events), "
               "not sepsis-vs-healthy."),
    ),
    "GSE69683": dict(
        label="Asthma, U-BIOPRED (peripheral blood)",
        platform="GPL13158",              # Affymetrix HT HG-U133+ PM
        expected_scale="log_normalized",
        outcome_hint=r"cohort|severity|group|status",
        notes=("87 healthy + 77 mild-moderate + 246 severe. Contrast: severe "
               "vs mild-moderate (published AUC around 0.70)."),
    ),
    "GSE13355": dict(
        label="Psoriasis, lesional vs uninvolved skin",
        platform="GPL570",                # Affymetrix U133 Plus 2.0
        expected_scale="log_normalized",
        outcome_hint=r"type|status|tissue|group",
        notes=("180 samples: NN normal control, PN uninvolved, PP involved. "
               "PAIRED (PN/PP from the same patient) -- CV folds MUST split on "
               "patient or the paired samples leak across the split. Known "
               "batch effect. Use as a positive control, not a discriminator."),
    ),
    "GSE107994": dict(
        label="Tuberculosis, Leicester (whole blood RNA-seq)",
        platform=None,                    # RNA-seq: symbols usually in suppl
        expected_scale="counts",
        outcome_hint=r"group|disease|state|status",
        notes=("53 active TB + 72 latent + 50 controls. Contrast: LTBI vs "
               "active. The progressor subset is far too small (~11) to use."),
    ),
    "GSE20194": dict(
        label="Breast cancer, MAQC-II neoadjuvant chemotherapy",
        platform="GPL96",                 # Affymetrix HG-U133A
        expected_scale="log_normalized",
        outcome_hint=r"pcr|response|residual|rcb|outcome",
        notes=("280 samples. The MAQC-II reference set, assembled specifically "
               "to benchmark prediction methods -- pathological complete "
               "response (pCR) vs residual disease after neoadjuvant "
               "chemotherapy. A genuinely hard contrast with published AUCs "
               "well short of ceiling."),
    ),
    "GSE25066": dict(
        label="Breast cancer, neoadjuvant taxane-anthracycline (pCR)",
        platform="GPL96",
        expected_scale="log_normalized",
        outcome_hint=r"pcr|response|residual|rcb|outcome",
        notes=("510 samples. Independent, larger replication of the same pCR "
               "vs residual-disease contrast as GSE20194."),
    ),
    "GSE16879": dict(
        label="IBD, infliximab response (mucosal biopsy)",
        platform="GPL570",
        expected_scale="log_normalized",
        outcome_hint=r"response|responder|infliximab|remission",
        notes=("135 samples. Predicting anti-TNF response from pre-treatment "
               "mucosa -- a problem that has resisted prediction for years, so "
               "the signal should be weak. Closest structural match to "
               "IMvigor210: does this patient respond to this biologic."),
    ),
    "GSE91061": dict(
        label="Melanoma, nivolumab response (RNA-seq)",
        platform=None,
        expected_scale="counts",
        outcome_hint=r"response|benefit|recist|therapy",
        notes=("111 samples. Checkpoint-blockade response, the same problem "
               "class as IMvigor210 in a different tumour type. RNA-seq, so "
               "expression is likely in a supplementary file."),
    ),
    "GSE59867": dict(
        label="Myocardial infarction, outcome (whole blood)",
        platform="GPL6244",
        expected_scale="log_normalized",
        outcome_hint=r"outcome|status|group|diagnosis|time",
        notes=("438 samples. Large cardiac cohort; prognosis from blood is a "
               "weak-signal problem. Contrast must be fixed from the metadata "
               "before running."),
    ),
    "GSE57945": dict(
        label="Crohn's disease, RISK cohort (ileal RNA-seq)",
        platform=None,
        expected_scale="counts",          # UNVERIFIED -- check the detection
        outcome_hint=r"diagnosis|group|disease",
        notes=("322 samples; CD 204 vs non-IBD control 42, plus UC. The authors "
               "issued a correction excluding 26 CD, 12 UC and 2 control "
               "samples -- apply it. Whether this ships counts or normalized "
               "values is UNVERIFIED; trust the detected scale."),
    ),
    "GSE92415": dict(
        label="Ulcerative colitis, golimumab response (baseline mucosa)",
        platform="GPL13158",             # Affymetrix HT HG-U133+ PM
        expected_scale="log_normalized",
        outcome_hint=r"wk6response|response|treatment|visit",
        notes=("Use baseline golimumab-treated biopsies only. GEO metadata "
               "give 32 week-6 responders and 27 non-responders. Patient IDs "
               "are in the 'subject' characteristic."),
    ),
    "GSE206285": dict(
        label="Ulcerative colitis, ustekinumab response (baseline mucosa)",
        platform="GPL13158",             # Affymetrix HT HG-U133+ PM
        expected_scale="log_normalized",
        outcome_hint=r"clinical_remission|mucosal_healing|treatment|visit",
        notes=("Use the 364 baseline ustekinumab-treated biopsies. The fixed "
               "endpoint is clinical remission at week 8: 49 remission and "
               "315 no remission. Keep donor ID and treatment dose."),
    ),
    "GSE19442": dict(
        label="Tuberculosis, South Africa external holdout (whole blood)",
        platform="GPL6947",              # Illumina HumanHT-12 V3.0
        expected_scale="linear_background_corrected",
        outcome_hint=r"disease|group|status|title",
        notes=("External holdout for panels locked in GSE107994. All 51 "
               "participants were sampled before treatment: 20 active PTB "
               "and 31 latent TB. Do not use it to rank methods by nested CV."),
    ),
    "GSE101794": dict(
        label="Paediatric Crohn disease, treatment-naive ileum",
        platform=None,                    # RNA-seq TPM in per-sample files
        expected_scale="log_normalized",
        outcome_hint=r"diagnosis|disease|group|title",
        notes=("304 treatment-naive ileal samples: 254 Crohn disease and 50 "
               "non-IBD. GEO supplies one TPM file per sample. Build "
               "expression.tsv with benchmarks/gse101794_build_expression.R; "
               "do not apply TMM to TPM."),
    ),
}


# =============================================================================
#  Download helpers
# =============================================================================

def geo_series_stem(accession):
    """
    GEO nests series by a truncated accession: GSE65682 lives under GSE65nnn.
    The rule is 'replace the last three digits with nnn', with short accessions
    becoming GSEnnn.
    """
    digits = accession[3:]
    return "GSE" + (digits[:-3] + "nnn" if len(digits) > 3 else "nnn")


def geo_platform_stem(platform):
    """Same nesting rule for platforms: GPL13667 -> GPL13nnn."""
    digits = platform[3:]
    return "GPL" + (digits[:-3] + "nnn" if len(digits) > 3 else "nnn")


def series_matrix_url(accession):
    return (f"{GEO_BASE}/series/{geo_series_stem(accession)}/{accession}"
            f"/matrix/{accession}_series_matrix.txt.gz")


def supplementary_dir_url(accession):
    return f"{GEO_BASE}/series/{geo_series_stem(accession)}/{accession}/suppl/"


def platform_annotation_url(platform):
    return (f"{GEO_BASE}/platforms/{geo_platform_stem(platform)}/{platform}"
            f"/annot/{platform}.annot.gz")


def download(url, destination, force=False):
    """
    Fetch `url` to `destination`, skipping if it already exists. Returns the
    path, or None if the resource is absent (404) -- absence is an expected
    outcome for several of these URLs (RNA-seq series often have no expression
    table; not every platform has an .annot.gz), so it is not an error.
    """
    destination = Path(destination)
    if destination.exists() and not force:
        size_mb = destination.stat().st_size / 1024 ** 2
        print(f"    [cached] {destination.name} ({size_mb:.1f} MB)")
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    print(f"    downloading {url}")

    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=120) as response, \
             open(destination, "wb") as handle:
            shutil.copyfileobj(response, handle)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            print(f"    [absent] {url} (404)")
            return None
        raise
    except Exception:
        # Never leave a truncated file behind: a partial .gz would be cached and
        # then fail to parse on the next run with a confusing error.
        if destination.exists():
            destination.unlink()
        raise

    size_mb = destination.stat().st_size / 1024 ** 2
    print(f"    saved {destination.name} ({size_mb:.1f} MB)")
    return destination


def list_supplementary_files(accession):
    """Scrape the suppl/ directory index for filenames."""
    url = supplementary_dir_url(accession)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            html = response.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError:
        return []
    names = re.findall(r'href="([^"?/][^"]*)"', html)
    return sorted({n for n in names if not n.startswith("..")})


# =============================================================================
#  Series matrix parsing
# =============================================================================

def open_maybe_gzip(path):
    """Text handle for a file that may or may not be gzipped."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_series_matrix(path):
    """
    Split a GEO series matrix into (expression, metadata).

    The file is two documents concatenated: '!'-prefixed header lines, then a
    tab-delimited table between the begin/end markers. RNA-seq series often
    carry the header but an EMPTY table, in which case expression is None and
    the counts have to come from the supplementary files instead.
    """
    sample_lines = {}
    table_start_line = None

    with open_maybe_gzip(path) as handle:
        for line_number, line in enumerate(handle):
            if line.startswith("!series_matrix_table_begin"):
                table_start_line = line_number
                break
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key[1:]                       # strip the leading '!'
            values = [v.strip().strip('"') for v in rest.split("\t")]
            # Characteristics keys repeat, one line per attribute; number them.
            sample_lines.setdefault(key, []).append(values)

    if table_start_line is None:
        raise ValueError(f"No series_matrix_table_begin marker in {path}")

    # --- Metadata ------------------------------------------------------------
    accessions = sample_lines.get("Sample_geo_accession", [[]])[0]
    metadata = pd.DataFrame(index=accessions)
    metadata.index.name = "sample"

    for key, blocks in sample_lines.items():
        if key == "Sample_geo_accession":
            continue
        for block_index, values in enumerate(blocks, start=1):
            if len(values) != len(accessions):
                continue
            suffix = "" if len(blocks) == 1 else f"_{block_index}"
            metadata[f"{key}{suffix}"] = values

    # Characteristics arrive as "field: value"; split them into real columns.
    # This is where the outcome variable lives, so it needs to be usable.
    for column in [c for c in metadata.columns
                   if c.startswith("Sample_characteristics")]:
        entries = metadata[column].astype(str)
        keys = entries.str.extract(r"^\s*([^:]+?)\s*:", expand=False)
        if keys.notna().mean() < 0.8:
            continue                        # not a key: value column
        field = keys.dropna().mode()
        if field.empty:
            continue
        field_name = re.sub(r"\W+", "_", field.iloc[0].strip().lower())
        metadata[f"ch_{field_name}"] = (
            entries.str.replace(r"^\s*[^:]+?\s*:\s*", "", regex=True).str.strip()
        )

    # --- Expression table ----------------------------------------------------
    expression = pd.read_csv(
        path, sep="\t", skiprows=table_start_line + 1, index_col=0,
        comment="!", low_memory=False,
        compression="gzip" if str(path).endswith(".gz") else None,
    )
    expression = expression.dropna(how="all")
    expression = expression.apply(pd.to_numeric, errors="coerce")

    if expression.shape[0] == 0 or expression.shape[1] == 0:
        expression = None

    return expression, metadata


# =============================================================================
#  Scale detection
# =============================================================================

def detect_scale(expression, sample_columns=200):
    """
    Decide whether a matrix is raw counts, log-scale, or linear-normalized, by
    looking at the values rather than at the filename.

    Heuristics, in order:
      * all values integral and the maximum is large   -> counts
      * negative values and maximum >= 100             -> background-corrected
                                                           linear array signal
      * remaining maximum < 30                         -> log_normalized
      * otherwise                                      -> linear_normalized
    """
    sample = expression.iloc[:, :min(sample_columns, expression.shape[1])]
    values = sample.to_numpy(dtype="float64", copy=False)
    finite = values[pd.notna(values)]
    if finite.size == 0:
        return "unknown", {}

    maximum = float(finite.max())
    minimum = float(finite.min())
    integral = bool((finite == finite.round()).all())
    negatives = bool((finite < 0).any())

    if integral and maximum > 100:
        scale = "counts"
    elif negatives and maximum >= 100:
        scale = "linear_background_corrected"
    elif maximum < 30:
        scale = "log_normalized"
    else:
        scale = "linear_normalized"

    return scale, dict(max=maximum, min=minimum, all_integers=integral,
                       has_negatives=negatives)


# =============================================================================
#  Probe -> symbol mapping
# =============================================================================

def load_platform_symbols(annotation_path):
    """
    Extract probe -> gene symbol from a GPL .annot.gz. The header is a block of
    '#'-prefixed comment lines followed by a tab-delimited table whose symbol
    column is named 'Gene symbol' (spelling varies by platform vintage).
    """
    with open_maybe_gzip(annotation_path) as handle:
        lines = handle.readlines()

    header_index = next((i for i, l in enumerate(lines)
                         if not l.startswith(("#", "^", "!"))), None)
    if header_index is None:
        return {}

    table = pd.read_csv(io.StringIO("".join(lines[header_index:])),
                        sep="\t", low_memory=False)

    id_column = table.columns[0]
    symbol_column = next(
        (c for c in table.columns
         if re.fullmatch(r"gene[\s_]*symbol", c.strip(), flags=re.I)), None)
    if symbol_column is None:
        return {}

    mapping = (table[[id_column, symbol_column]]
               .dropna()
               .astype(str)
               .set_index(id_column)[symbol_column]
               .to_dict())
    # Multi-mapping probes ("A///B") are ambiguous; drop rather than pick one.
    return {probe: symbol for probe, symbol in mapping.items()
            if symbol and symbol != "---" and "///" not in symbol}


def collapse_to_symbols(expression, probe_to_symbol):
    """
    Map probes to symbols and resolve duplicates by keeping the probe with the
    highest mean signal.

    Averaging probes for the same gene is the tempting alternative but it is
    worse: platforms carry non-responsive probes whose flat signal drags the
    average toward background, damping real differences. Taking the strongest
    probe is the conventional choice.
    """
    symbols = pd.Series(expression.index.map(probe_to_symbol),
                        index=expression.index)
    keep = symbols.notna()
    dropped = int((~keep).sum())

    mapped = expression.loc[keep].copy()
    mapped["__symbol__"] = symbols[keep].values
    mapped["__mean__"] = expression.loc[keep].mean(axis=1).values

    collapsed = (mapped
                 .sort_values("__mean__", ascending=False)
                 .drop_duplicates("__symbol__", keep="first")
                 .set_index("__symbol__")
                 .drop(columns="__mean__"))
    collapsed.index.name = "gene"
    return collapsed, dropped


# =============================================================================
#  Per-dataset driver
# =============================================================================

def fetch_dataset(accession, outdir, force=False, inspect_only=False):
    config = DATASETS[accession]
    target = Path(outdir) / accession
    raw = target / "raw"

    print(f"\n=== {accession}: {config['label']} ===")
    print(f"  {config['notes']}")

    matrix_path = download(series_matrix_url(accession),
                           raw / f"{accession}_series_matrix.txt.gz",
                           force=force)
    if matrix_path is None:
        print("  !! series matrix unavailable; cannot continue")
        return None

    print("  parsing series matrix...")
    expression, metadata = parse_series_matrix(matrix_path)

    # --- Report the characteristics so the outcome column is chosen by eye ---
    characteristic_columns = [c for c in metadata.columns if c.startswith("ch_")]
    print(f"  {metadata.shape[0]} samples | "
          f"{len(characteristic_columns)} parsed characteristics:")
    for column in characteristic_columns:
        levels = metadata[column].value_counts()
        preview = ", ".join(f"{k}={v}" for k, v in levels.head(4).items())
        flag = " <-- candidate outcome" if re.search(
            config["outcome_hint"], column, flags=re.I) else ""
        print(f"      {column:36s} {len(levels):4d} levels  [{preview}]{flag}")

    if inspect_only:
        return None

    # --- Expression ---------------------------------------------------------
    manifest = dict(accession=accession, label=config["label"],
                    n_samples=int(metadata.shape[0]),
                    platform=config["platform"],
                    expected_scale=config["expected_scale"],
                    notes=config["notes"])

    if expression is None:
        # Typical for RNA-seq series: the matrix carries metadata only.
        print("  series matrix has NO expression table (usual for RNA-seq).")
        files = list_supplementary_files(accession)
        print(f"  {len(files)} supplementary file(s) on the server:")
        for name in files:
            print(f"      {name}")
        manifest.update(scale="unknown", n_genes=0,
                        supplementary_files=files,
                        action_required=(
                            "Expression is in a supplementary file. Download "
                            "the counts matrix, then re-run with "
                            "--from-supplementary <filename>."))
    else:
        scale, statistics = detect_scale(expression)
        print(f"  expression: {expression.shape[0]} rows x "
              f"{expression.shape[1]} samples")
        print(f"  detected scale: {scale}  "
              f"(max={statistics['max']:.3g}, "
              f"integers={statistics['all_integers']}, "
              f"negatives={statistics['has_negatives']})")

        if scale != config["expected_scale"]:
            print(f"  !! SCALE MISMATCH: expected {config['expected_scale']}, "
                  f"detected {scale}. Trust the detection and check the GEO "
                  f"record before proceeding.")

        # --- Probe -> symbol ------------------------------------------------
        dropped = None
        if config["platform"]:
            annotation = download(
                platform_annotation_url(config["platform"]),
                raw / f"{config['platform']}.annot.gz", force=force)
            if annotation is None:
                print(f"  !! no .annot.gz for {config['platform']}; rows remain "
                      f"PROBE IDs. The biology pillar needs symbols -- map them "
                      f"before running the benchmark.")
            else:
                probe_to_symbol = load_platform_symbols(annotation)
                if not probe_to_symbol:
                    print("  !! annotation had no usable 'Gene symbol' column; "
                          "rows remain probe IDs.")
                else:
                    before = expression.shape[0]
                    expression, dropped = collapse_to_symbols(
                        expression, probe_to_symbol)
                    print(f"  mapped probes to symbols: {before} probes -> "
                          f"{expression.shape[0]} genes "
                          f"({dropped} probes unmapped)")

        manifest.update(scale=scale, n_genes=int(expression.shape[0]),
                        value_stats=statistics,
                        unmapped_probes=dropped,
                        row_identifier=("gene_symbol" if config["platform"]
                                        and dropped is not None else "probe_id"))

        target.mkdir(parents=True, exist_ok=True)
        expression_path = target / "expression.tsv"
        expression.to_csv(expression_path, sep="\t")
        print(f"  wrote {expression_path}")

    target.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(target / "metadata.tsv", sep="\t")
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"  wrote {target / 'metadata.tsv'}")
    print(f"  wrote {target / 'manifest.json'}")
    return manifest


def main():
    parser = argparse.ArgumentParser(
        description="Fetch GEO datasets for the GeneSelectR benchmark.")
    parser.add_argument("accessions", nargs="*", help="GSE accessions")
    parser.add_argument("--all", action="store_true", help="fetch every dataset")
    parser.add_argument("--list", action="store_true", help="list the registry")
    parser.add_argument("--outdir", default="data", help="output directory")
    parser.add_argument("--force", action="store_true", help="re-download")
    parser.add_argument("--inspect", action="store_true",
                        help="print characteristics only; write nothing")
    arguments = parser.parse_args()

    if arguments.list:
        print(f"{'accession':12s} {'expected scale':18s} label")
        for accession, config in DATASETS.items():
            print(f"{accession:12s} {config['expected_scale']:18s} "
                  f"{config['label']}")
        return 0

    targets = list(DATASETS) if arguments.all else arguments.accessions
    if not targets:
        parser.print_help()
        return 1

    unknown = [a for a in targets if a not in DATASETS]
    if unknown:
        return f"Not in the registry: {', '.join(unknown)}"

    manifests = []
    for accession in targets:
        try:
            result = fetch_dataset(accession, arguments.outdir,
                                   force=arguments.force,
                                   inspect_only=arguments.inspect)
            if result:
                manifests.append(result)
        except Exception as error:
            print(f"  !! {accession} FAILED: {type(error).__name__}: {error}")

    if manifests:
        print(f"\n{'=' * 70}\nSummary\n{'=' * 70}")
        print(f"{'accession':12s} {'samples':>8s} {'genes':>8s} {'scale':20s}")
        for manifest in manifests:
            print(f"{manifest['accession']:12s} {manifest['n_samples']:>8d} "
                  f"{manifest['n_genes']:>8d} {manifest['scale']:20s}")
        summary = Path(arguments.outdir) / "manifest_all.json"
        summary.write_text(json.dumps(manifests, indent=2))
        print(f"\nwrote {summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
