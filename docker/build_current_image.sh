#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

IMAGE_REPOSITORY="${UNISCFLOW_IMAGE_REPOSITORY:-uniscflow}"
TARGET_PLATFORM="${UNISCFLOW_DOCKER_PLATFORM:-linux/amd64}"

if ! command -v git >/dev/null 2>&1; then
  echo "ERROR: git is required to bind the image to an exact source revision." >&2
  exit 1
fi
if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: docker is not available on PATH." >&2
  exit 1
fi
if [[ -n "$(git status --porcelain)" ]]; then
  echo "ERROR: refusing to build a provenance-labelled image from a dirty checkout." >&2
  echo "Commit or remove local changes, then run this script again." >&2
  exit 1
fi

GIT_COMMIT="$(git rev-parse HEAD)"
SHORT_COMMIT="$(git rev-parse --short=7 HEAD)"
VERSION="$(sed -n 's/^__version__ = "\([^"]*\)"/\1/p' uniscflow.py | head -n 1)"
BUILD_DATE="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"

if [[ -z "${VERSION}" ]]; then
  echo "ERROR: could not read __version__ from uniscflow.py." >&2
  exit 1
fi

LATEST_TAG="${IMAGE_REPOSITORY}:latest"
VERSION_TAG="${IMAGE_REPOSITORY}:${VERSION}"
COMMIT_TAG="${IMAGE_REPOSITORY}:${SHORT_COMMIT}"

docker build --pull --platform "${TARGET_PLATFORM}" \
  --build-arg "UNISCFLOW_VERSION=${VERSION}" \
  --build-arg "UNISCFLOW_GIT_COMMIT=${GIT_COMMIT}" \
  --build-arg "UNISCFLOW_BUILD_DATE=${BUILD_DATE}" \
  --tag "${LATEST_TAG}" \
  --tag "${VERSION_TAG}" \
  --tag "${COMMIT_TAG}" \
  .

bash docker/verify_image.sh "${COMMIT_TAG}" "${GIT_COMMIT}" "${VERSION}"

printf 'Built and verified:\n  %s\n  %s\n  %s\nRevision: %s\n' \
  "${LATEST_TAG}" "${VERSION_TAG}" "${COMMIT_TAG}" "${GIT_COMMIT}"
