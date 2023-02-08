# Nextstrain base with procps

FROM nextstrain/base:build-20230208T000021Z

RUN apt-get update && apt-get install -y procps

ENTRYPOINT ["augur"]
