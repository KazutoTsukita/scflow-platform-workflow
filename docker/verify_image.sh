#!/usr/bin/env bash
set -euo pipefail

IMAGE="${1:-uniscflow:latest}"
EXPECTED_REVISION="${2:-}"
EXPECTED_VERSION="${3:-}"
TARGET_PLATFORM="${UNISCFLOW_DOCKER_PLATFORM:-linux/amd64}"

HELP_OUTPUT="$(mktemp)"
PLAN_OUTPUT="$(mktemp)"
trap 'rm -f "${HELP_OUTPUT}" "${PLAN_OUTPUT}"' EXIT

docker run --rm --platform "${TARGET_PLATFORM}" "${IMAGE}" --version
docker run --rm --platform "${TARGET_PLATFORM}" "${IMAGE}" > "${HELP_OUTPUT}"
grep -q '^usage: uniscflow' "${HELP_OUTPUT}"
docker run --rm --platform "${TARGET_PLATFORM}" "${IMAGE}" --mode plan > "${PLAN_OUTPUT}"
grep -q '^Action: plan$' "${PLAN_OUTPUT}"

docker run --rm --platform "${TARGET_PLATFORM}" --entrypoint bash "${IMAGE}" -c '
  set -euo pipefail
  for executable in STAR featureCounts samtools fasterq-dump pigz wget parallel Rscript salmon python; do
    command -v "${executable}" >/dev/null
  done
  test -s /workflow/UNISCFLOW_IMAGE_PROVENANCE
  test -d /workflow/profiles/platforms
  test "$(find /workflow/profiles/platforms -maxdepth 1 -name "*.json" -type f | wc -l | tr -d " ")" = 25
  test "$(id -u)" != 0
'

IMAGE_REVISION="$(docker image inspect \
  --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}' \
  "${IMAGE}")"
IMAGE_VERSION="$(docker image inspect \
  --format '{{ index .Config.Labels "org.opencontainers.image.version" }}' \
  "${IMAGE}")"
IMAGE_USER="$(docker image inspect --format '{{ .Config.User }}' "${IMAGE}")"

if [[ -n "${EXPECTED_REVISION}" && "${IMAGE_REVISION}" != "${EXPECTED_REVISION}" ]]; then
  echo "ERROR: image revision ${IMAGE_REVISION} does not match ${EXPECTED_REVISION}." >&2
  exit 1
fi
if [[ -n "${EXPECTED_VERSION}" && "${IMAGE_VERSION}" != "${EXPECTED_VERSION}" ]]; then
  echo "ERROR: image version ${IMAGE_VERSION} does not match ${EXPECTED_VERSION}." >&2
  exit 1
fi
if [[ -z "${IMAGE_USER}" || "${IMAGE_USER}" == "root" || "${IMAGE_USER}" == "0" ]]; then
  echo "ERROR: image is not configured to run as a non-root user." >&2
  exit 1
fi

printf 'Verified image: %s\nVersion: %s\nRevision: %s\nUser: %s\n' \
  "${IMAGE}" "${IMAGE_VERSION}" "${IMAGE_REVISION}" "${IMAGE_USER}"
