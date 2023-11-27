FROM alpine:3

COPY target/x86_64-unknown-linux-musl/release/rebar /usr/bin/rebar

RUN mkdir -p /rebar/dataset

ADD ./data /rebar/data
COPY ./dataset/toy* /rebar/dataset

WORKDIR "/rebar"

RUN set -eux
