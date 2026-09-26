FROM ubuntu:22.04

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
    jq \
    less \
    procps \
    python3 \
    sed \
    tar \
    wget \
    xz-utils \
    && rm -rf /var/lib/apt/lists/*

COPY inspect_cellranger_chemistry.py /usr/local/bin/inspect_cellranger_chemistry
RUN chmod +x /usr/local/bin/inspect_cellranger_chemistry

WORKDIR /opt
CMD ["/bin/bash"]
