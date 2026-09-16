# syntax=docker/dockerfile:1

# Debian slim base: apt ships hmmer for both amd64 and arm64, so this Dockerfile
# builds natively on Apple Silicon (via `docker buildx build --platform linux/arm64`)
# and on x86 Ubuntu/CI, with no per-arch branching. This sidesteps the bioconda
# `anarci` package, which only publishes an osx-64 build and forces the Rosetta
# workaround documented in README.md's addendum.
FROM python:3.10-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    hmmer \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/oxpig/ANARCI.git /opt/ANARCI \
    && cd /opt/ANARCI \
    && python3 setup.py install
ENV PATH="${PATH}:/opt/ANARCI/bin"

# `setup.py install` shells out to easy_install to resolve AbNumber's declared
# deps (pandas, anarcii), which chokes on modern sdist-less pandas releases.
# --no-deps installs just the AbNumber package itself; requirements.txt below
# covers pandas/biopython, and anarcii is unrelated to our ANARCI (oxpig) install.
RUN git clone --depth 1 https://github.com/prihoda/AbNumber.git /opt/AbNumber \
    && pip install --no-cache-dir --no-deps /opt/AbNumber

WORKDIR /app

COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

COPY pyproject.toml uv.lock MANIFEST.in README.md ./
COPY src/ src/

RUN uv pip install --system --no-cache .

# Test-only deps (the `dev` extra — not needed to *use* anarci-toolz, only to
# run its test suite, which is what this image's CI job does via
# `--entrypoint pytest`).
RUN uv pip install --system --no-cache .[dev]

COPY pytest.ini ./
COPY test_files/ test_files/
COPY tests/ tests/

ENTRYPOINT ["anarci-toolz"]
