# FlashFry image: FlashFry ships baked in, so end users never install or
# download it themselves. Pinned to 1.15 to match the index-build command
# documented in assets/build_reference_index.md.

FROM openjdk:8-jdk-slim

ARG FLASHFRY_VERSION=1.15
RUN apt-get update && apt-get install -y --no-install-recommends wget && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/flashfry && \
    wget -O /opt/flashfry/FlashFry-assembly-${FLASHFRY_VERSION}.jar \
    https://github.com/mckennalab/FlashFry/releases/download/${FLASHFRY_VERSION}/FlashFry-assembly-${FLASHFRY_VERSION}.jar && \
    ln -s /opt/flashfry/FlashFry-assembly-${FLASHFRY_VERSION}.jar /opt/flashfry/FlashFry-assembly-1.15.jar

WORKDIR /work
