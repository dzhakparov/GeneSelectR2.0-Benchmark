# GeneSelectR analysis code

This branch contains the scripts used for GeneSelectR method development,
configuration comparisons, ablation studies, validation, biological
assessment and figure generation. Generated results and input data are not
committed.

Start with [analysis/README.md](analysis/README.md). The analysis order is
implemented in `analysis/run_all.sh`, and the common settings are stored in
`analysis/config.R`.

The package source is maintained on `codex/bioconductor-release`.
`analysis/bootstrap_package.sh` extracts package revision `621be0c1` into a
local, ignored `package/` directory before an analysis run. The
`GENESELECTR_PACKAGE_REF` environment variable can select another revision.
