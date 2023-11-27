FROM alpine:3

# install fontconfig
RUN apk --no-cache add fontconfig ttf-dejavu && \
    fc-cache -f -v

COPY target/x86_64-unknown-linux-musl/release/rebar /usr/bin/rebar

RUN mkdir -p /rebar/dataset

ADD ./data /rebar/data
COPY ./dataset/toy* /rebar/dataset

WORKDIR "/rebar"

RUN set -eux
