FROM mambaorg/micromamba:1.5.8

ARG UNISCFLOW_VERSION=unknown
ARG UNISCFLOW_GIT_COMMIT=unknown
ARG UNISCFLOW_BUILD_DATE=unknown

LABEL org.opencontainers.image.title="UniScFlow" \
      org.opencontainers.image.description="Public-accession to validated scRNA-seq mapping workflow" \
      org.opencontainers.image.source="https://github.com/KazutoTsukita/scflow-platform-workflow" \
      org.opencontainers.image.url="https://github.com/KazutoTsukita/scflow-platform-workflow" \
      org.opencontainers.image.documentation="https://github.com/KazutoTsukita/scflow-platform-workflow/blob/main/README.md" \
      org.opencontainers.image.authors="Kazuto Tsukita and contributors" \
      org.opencontainers.image.vendor="UniScFlow" \
      org.opencontainers.image.licenses="BSD-3-Clause" \
      org.opencontainers.image.version="${UNISCFLOW_VERSION}" \
      org.opencontainers.image.revision="${UNISCFLOW_GIT_COMMIT}" \
      org.opencontainers.image.created="${UNISCFLOW_BUILD_DATE}"

USER root
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    bzip2 \
    ca-certificates \
    coreutils \
    curl \
    findutils \
    gawk \
    gzip \
    git \
    jq \
    procps \
    sed \
    tar \
    tini \
    util-linux \
    xz-utils \
    && rm -rf /var/lib/apt/lists/*

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.10 \
    tomli \
    pip \
    pandas \
    r-base \
    r-readr \
    r-dplyr \
    r-stringr \
    parallel \
    pigz \
    wget \
    sra-tools \
    star=2.7.10b \
    samtools \
    subread \
    salmon \
    && micromamba clean -a -y

ENV PATH=/workflow/bin:/opt/conda/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

WORKDIR /workflow
COPY . /workflow
RUN printf 'version=%s\nrevision=%s\nbuild_date=%s\n' \
      "${UNISCFLOW_VERSION}" "${UNISCFLOW_GIT_COMMIT}" "${UNISCFLOW_BUILD_DATE}" \
      > /workflow/UNISCFLOW_IMAGE_PROVENANCE \
    && python -m pip install --no-deps --no-cache-dir /workflow \
    && mkdir -p /data \
    && chown -R "${MAMBA_USER}:${MAMBA_USER}" /workflow /data
USER ${MAMBA_USER}
WORKDIR /data

ENTRYPOINT ["tini", "--", "/workflow/bin/uniscflow"]
CMD ["--help"]
