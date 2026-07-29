#!/usr/bin/env bash
#
# Run the workflow twice on the same input -- once with the classic `pairtools`
# backend, once with `pairtools_parquet` -- into two isolated working directories,
# so the per-rule `benchmark:` TSVs of the two arms can be compared.
#
#   benchmarking/run_benchmark.sh [-c CORES] [-f CONFIGFILE] [-o OUTDIR] [-a ARM] [-- ...]
#
# Anything after `--` is passed straight through to snakemake.
#
set -eo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

CORES=4
CONFIGFILE="${REPO_ROOT}/config/config.yml"
OUTDIR="${REPO_ROOT}/bench"
ARMS="pairtools parquet"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--cores)      CORES="$2"; shift 2 ;;
        -f|--configfile) CONFIGFILE="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"; shift 2 ;;
        -o|--outdir)     OUTDIR="$2"; shift 2 ;;
        -a|--arm)        ARMS="$2"; shift 2 ;;
        --)              shift; break ;;
        -h|--help)
            sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
EXTRA_ARGS=("$@")

mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"

# GNU time gives a total including snakemake's own scheduling overhead, which the
# per-rule benchmark TSVs do not capture. It is not always installed (the bash
# builtin `time` cannot write to a file), so fall back to timing the run ourselves.
GNU_TIME="$(command -v gtime || command -v /usr/bin/time || true)"
if [[ -z "${GNU_TIME}" ]]; then
    echo "note: GNU time not found; total wall clock will be measured with \$SECONDS" >&2
fi

run_timed() {  # run_timed <outfile> <cmd...>
    local outfile="$1"; shift
    if [[ -n "${GNU_TIME}" ]]; then
        "${GNU_TIME}" -v -o "${outfile}" "$@"
    else
        local start=${SECONDS}
        "$@"
        local status=$?
        local elapsed=$((SECONDS - start))
        printf '\tElapsed (wall clock) time (h:mm:ss or m:ss): %d:%02d:%02d\n' \
            $((elapsed / 3600)) $(((elapsed % 3600) / 60)) $((elapsed % 60)) > "${outfile}"
        return ${status}
    fi
}

# Inputs referenced by the config are relative to the repo root, so link them into
# each arm's working directory rather than copying. The bwa index is built once and
# shared between the arms through snakemake's between-workflow cache (rule bwaindex
# is marked `cache: True`).
export SNAKEMAKE_OUTPUT_CACHE="${OUTDIR}/snakemake-cache"
mkdir -p "${SNAKEMAKE_OUTPUT_CACHE}"

for arm in ${ARMS}; do
    workdir="${OUTDIR}/${arm}"
    mkdir -p "${workdir}"

    # Link every top-level entry that inputs may refer to (test data, custom genomes).
    for entry in "${REPO_ROOT}"/*; do
        name="$(basename "${entry}")"
        case "${name}" in
            bench|workflow|config|benchmarking|.git) continue ;;
        esac
        [[ -e "${workdir}/${name}" ]] || ln -s "${entry}" "${workdir}/${name}"
    done

    # A second configfile rather than --config: snakemake deep-merges configfiles,
    # so the rest of the `pairtools:` section (dedup_backend, memory, tmpdir) is kept,
    # whereas --config would replace the whole mapping.
    printf 'pairtools:\n    backend: %s\n' "${arm}" > "${workdir}/backend.yml"

    echo "=============================================================="
    echo "  arm: ${arm}   cores: ${CORES}   workdir: ${workdir}"
    echo "=============================================================="

    run_timed "${workdir}/total_time.txt" \
        snakemake \
            --snakefile "${REPO_ROOT}/workflow/Snakefile" \
            --configfile "${CONFIGFILE}" "${workdir}/backend.yml" \
            --directory "${workdir}" \
            --cores "${CORES}" \
            --use-conda \
            --cache \
            "${EXTRA_ARGS[@]}" \
        2>&1 | tee "${workdir}/snakemake.log"
done

echo
echo "Done. Compare with:"
echo "  python benchmarking/compare.py ${OUTDIR}"
echo "  python benchmarking/compare_outputs.py ${OUTDIR}"
