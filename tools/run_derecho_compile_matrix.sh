#!/usr/bin/env bash

# Compile the current committed source on Derecho with the documented NVHPC
# case profiles.  This is intentionally separate from the local readiness gate
# because it requires an authenticated `ssh derecho` session.
#
# Usage from the repository root:
#
#   tools/run_derecho_compile_matrix.sh
#
# Optional environment variables:
#
#   DERECHO_COMPILE_ROOT   Remote output directory.  Defaults to a timestamped
#                          directory under /glade/work/wchen/lesgo_versions.
#   BUILD_JOBS             Parallel build jobs passed to the case compile
#                          scripts.  Defaults to 8.
#   DERECHO_CASES          Space-separated case list.  Defaults to the four
#                          public presentation cases.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "${repo_root}"

if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "Working tree has uncommitted changes; commit first so Derecho tests a reproducible source snapshot." >&2
    exit 2
fi

commit="$(git rev-parse --short HEAD)"
timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
remote_root="${DERECHO_COMPILE_ROOT:-/glade/work/wchen/lesgo_versions/benchmarks/derecho_compile_matrix_${commit}_${timestamp}}"
build_jobs="${BUILD_JOBS:-8}"
cases_text="${DERECHO_CASES:-channel_flow adm_disk atm_line large_windfarm_3072x384x400_60turbines}"

read -r -a cases <<< "${cases_text}"
if [[ "${#cases[@]}" -eq 0 ]]; then
    echo "No cases requested." >&2
    exit 2
fi

ssh derecho bash -s <<REMOTE
set -euo pipefail
mkdir -p "${remote_root}"
rm -rf "${remote_root}/src"
mkdir -p "${remote_root}/src"
REMOTE

git archive HEAD | ssh derecho "tar -xf - -C '${remote_root}/src'"
ssh derecho "printf '%s\n' '${commit}' > '${remote_root}/src/LOCAL_COMMIT.txt'"

printf 'Derecho compile root: %s\n' "${remote_root}"
printf 'Commit: %s\n' "${commit}"

for case_name in "${cases[@]}"; do
    printf '\n== %s gpu ==\n' "${case_name}"
    ssh derecho bash -s -- "${remote_root}" "${case_name}" "${build_jobs}" <<'REMOTE'
set -euo pipefail
remote_root="$1"
case_name="$2"
build_jobs="$3"
src="${remote_root}/src"
log="${remote_root}/compile_${case_name}_gpu.log"

rm -rf "${src}/test-cases/${case_name}/build-derecho-gpu"
if (cd "${src}" && BUILD_JOBS="${build_jobs}" bash "test-cases/${case_name}/compile_derecho.sh" gpu > "${log}" 2>&1); then
    printf 'PASS %s %s\n' "${case_name}" "${log}"
    grep -E 'Built .*/lesgo-run-exe-gpu|Built target|Linking Fortran executable' "${log}" | tail -n 5 || true
else
    rc=$?
    printf 'FAIL %s rc=%s %s\n' "${case_name}" "${rc}" "${log}"
    tail -n 140 "${log}"
    exit "${rc}"
fi
REMOTE
done

printf '\nAll requested Derecho GPU compile profiles passed for %s.\n' "${commit}"
