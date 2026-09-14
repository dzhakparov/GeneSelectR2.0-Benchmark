#!/usr/bin/env bash
set -euo pipefail

repository_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$repository_root"

for ref in 621be0c1 db1f6528; do
  git rev-parse --verify --quiet "${ref}^{commit}" >/dev/null || {
    printf 'Missing package revision: %s\n' "$ref" >&2
    exit 1
  }
done

temporary_root="$(mktemp -d /private/tmp/geneselectr_package_revision.XXXXXX)"
trap 'rm -rf "$temporary_root"' EXIT
mkdir -p "$temporary_root/621" "$temporary_root/db1"
git archive 621be0c1 package/GeneSelectR/R | tar -x -C "$temporary_root/621"
git archive db1f6528 package/GeneSelectR/R | tar -x -C "$temporary_root/db1"

Rscript - "$temporary_root" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
required <- c(
    "gate_method", "subsample_scheme", "subsample_k_folds",
    "utility_method", "components", "score_formula", "calibration_mode",
    "calibration_n_permutations", "calibration_null_B", "use_cache",
    "n_cores"
)
for (ref in c("621", "db1")) {
    env <- new.env(parent = globalenv())
    files <- list.files(
        file.path(args[[1]], ref, "package", "GeneSelectR", "R"),
        full.names = TRUE
    )
    for (file in files) sys.source(file, env)
    available <- names(formals(env$geneselectr2_fit))
    missing <- setdiff(required, available)
    cat(sprintf("%s: missing reported-run arguments: %s\n", ref,
                if (length(missing)) paste(missing, collapse = ", ") else "none"))
    if (ref == "621" && length(missing)) quit(save = "no", status = 1L)
    if (ref == "db1" && !length(missing)) quit(save = "no", status = 1L)
}
RSCRIPT
