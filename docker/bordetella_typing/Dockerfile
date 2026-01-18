FROM ubuntu:24.04

WORKDIR /app

RUN apt-get update && \
    apt-get install -y \
        jq \
        curl \
        python3 && \
    apt-get clean

COPY ./dev/* /app/
COPY ./bin/bordetella_typing /bin/
RUN chmod 755 /bin/bordetella_typing

