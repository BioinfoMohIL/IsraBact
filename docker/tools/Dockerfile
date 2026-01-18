FROM ubuntu:22.04

RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends curl jq coreutils ca-certificates && \
    rm -rf /var/lib/apt/lists/*
