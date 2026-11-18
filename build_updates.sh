#!/usr/bin/env bash
# Build the updated images and push them to Docker Hub (evolbioinfo/<tool>:<version>).
#
# Usage:
#   ./build_updates.sh                 build and push every image listed in IMAGES
#   ./build_updates.sh samtools jq     only the listed tools
#   ./build_updates.sh --no-push ...   build only, do not push
#   ./build_updates.sh --dry-run ...   print the commands without running them
#
# Each image is built for linux/amd64 (the Dockerfiles download x86_64
# binaries). A failed build does not stop the others: failures are listed at
# the end, and each build's output is kept in build_logs/<tool>_<version>.log.
# Requires 'docker login' with an account that can push to evolbioinfo.

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
ORG="evolbioinfo"
PLATFORM="linux/amd64"
LOG_DIR="${REPO_ROOT}/build_logs"

# <tool directory>/<version directory>  (the version directory is the image tag)
IMAGES=(
    adapterremoval/v3.0.2
    alfred/v0.6.2
    bbmap/v40.05
    bcftools/v1.25
    checkm/v1.3.0
    fastqc/v0.13.0
    gubbins/v3.5
    hisat2/v2.2.3
    hyphy/v2.5.103
    igv/v2.19.8
    iqtree/v3.1.4
    irma/v1.4.0
    jq/v1.8.2
    kraken/v2.17.2
    label/v0.7.3
    lastal/v1654
    minimap2/v2.31
    mrbayes/v3.2.9
    multiqc/v1.35
    nanoplot/v1.48.0
    paml/v4.10.10
    pathphynder/v1.3.1
    phyml-sms/v1.8.3
    picard/v3.5.0
    predicthaplo/v2.1.4
    raxml-ng/v2.0.4
    samtools/v1.25
    schmutzi/v1.5.7
    seqkit/v2.14.0
    snakemake/v9.27.2
    spades/v4.3.0
    sphinx/v9.2.0
    sra-tools/v3.4.1
    sratoolkit/v3.4.1
    stringtie/v3.0.3
    treesimulator/v0.2.30
    treetime/v0.12.1
    trimgalore/v2.5.0
    wget/v1.25.0
)

PUSH=1
DRY_RUN=0
SELECTED=()
for arg in "$@"; do
    case "$arg" in
        --no-push) PUSH=0 ;;
        --dry-run) DRY_RUN=1 ;;
        -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
        -*) echo "Unknown option: $arg" >&2; exit 2 ;;
        *) SELECTED+=("$arg") ;;
    esac
done

# Keep only the requested tools, if any were given
TODO=()
if [ ${#SELECTED[@]} -eq 0 ]; then
    TODO=("${IMAGES[@]}")
else
    for tool in "${SELECTED[@]}"; do
        found=0
        for image in "${IMAGES[@]}"; do
            if [ "${image%%/*}" = "$tool" ]; then TODO+=("$image"); found=1; fi
        done
        [ $found -eq 1 ] || { echo "Unknown tool: $tool" >&2; exit 2; }
    done
fi

run() {
    if [ $DRY_RUN -eq 1 ]; then
        echo "+ $*"
    else
        "$@"
    fi
}

# Check everything before starting the (long) builds
for image in "${TODO[@]}"; do
    [ -f "${REPO_ROOT}/${image}/Dockerfile" ] || { echo "Missing ${image}/Dockerfile" >&2; exit 1; }
done
if [ $DRY_RUN -eq 0 ]; then
    docker info >/dev/null 2>&1 || { echo "Docker is not running" >&2; exit 1; }
    mkdir -p "$LOG_DIR"
fi

OK=()
FAILED=()
for image in "${TODO[@]}"; do
    tool="${image%%/*}"
    version="${image#*/}"
    tag="${ORG}/${tool}:${version}"
    log="${LOG_DIR}/${tool}_${version}.log"
    echo "=== ${tag}"

    if [ $DRY_RUN -eq 1 ]; then
        run docker build --platform "$PLATFORM" -t "$tag" "${REPO_ROOT}/${image}"
        [ $PUSH -eq 1 ] && run docker push "$tag"
        continue
    fi

    if ! docker build --platform "$PLATFORM" -t "$tag" "${REPO_ROOT}/${image}" >"$log" 2>&1; then
        echo "    build FAILED (see ${log})"
        FAILED+=("$tag (build)")
        continue
    fi
    echo "    built"

    if [ $PUSH -eq 1 ]; then
        if ! docker push "$tag" >>"$log" 2>&1; then
            echo "    push FAILED (see ${log})"
            FAILED+=("$tag (push)")
            continue
        fi
        echo "    pushed"
    fi
    OK+=("$tag")
done

[ $DRY_RUN -eq 1 ] && exit 0

echo
echo "Succeeded: ${#OK[@]}"
for tag in "${OK[@]}"; do echo "  $tag"; done
if [ ${#FAILED[@]} -gt 0 ]; then
    echo "Failed: ${#FAILED[@]}"
    for tag in "${FAILED[@]}"; do echo "  $tag"; done
    exit 1
fi
