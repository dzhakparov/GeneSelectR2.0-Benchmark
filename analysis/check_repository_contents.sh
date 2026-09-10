#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

allowed='(^|/)(LICENSE|README[.]md|[.]gitignore|[.]gitattributes)$|[.](R|py|sh|md|png|jpg|jpeg|svg)$'
unexpected=()

while IFS= read -r path; do
  if [[ ! "$path" =~ $allowed ]]; then
    unexpected+=("$path")
  fi
done < <(git ls-files)

if (( ${#unexpected[@]} > 0 )); then
  printf 'Tracked files outside the code, documentation and image policy:\n' >&2
  printf '  %s\n' "${unexpected[@]}" >&2
  exit 1
fi

printf 'Tracked-file policy passed. No data or generated result files are tracked.\n'
