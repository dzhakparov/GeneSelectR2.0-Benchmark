#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

package_ref="${GENESELECTR_PACKAGE_REF:-621be0c1}"
package_dir="package/GeneSelectR"

if [[ -d "$package_dir/R" ]]; then
  printf 'Package source already present: %s\n' "$package_dir"
  exit 0
fi

if ! git rev-parse --verify --quiet "${package_ref}^{commit}" >/dev/null; then
  printf 'Package ref is unavailable: %s\n' "$package_ref" >&2
  printf 'Fetch the package branch or set GENESELECTR_PACKAGE_REF.\n' >&2
  exit 1
fi

git archive "$package_ref" package/GeneSelectR | tar -x
printf 'Package source extracted from %s\n' "$package_ref"
