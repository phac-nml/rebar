FROM alpine:3

COPY target/x86_64-unknown-linux-musl/release/rebar /usr/bin/rebar

WORKDIR "/rebar"

RUN set -eux
