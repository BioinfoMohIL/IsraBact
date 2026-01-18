FROM ubuntu:24.04

WORKDIR /app

RUN apt-get update && \
    apt-get install -y \
        jq \
        curl \
        python3 python3-pandas && \
    apt-get clean

COPY ./dev/* /app/
COPY ./bin/neisseria_typing /bin/
RUN chmod 755 /bin/neisseria_typing

# ENTRYPOINT ["python3", "/bin/"]
