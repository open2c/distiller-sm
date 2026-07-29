#!/usr/bin/env bash
#
# Snakemake runs this inside the freshly created environment.
#
# pairtools ships a Cython extension (pairtools/lib/parse_pysam) that links against
# pysam. Under pip's default build isolation it is compiled against a throwaway copy
# of pysam in a temporary overlay, and the resulting .so keeps pointing at that
# directory -- so every `pairtools parse` then dies with
#
#   ImportError: /tmp/pip-build-env-*/overlay/.../pysam/libchtslib...so:
#                cannot open shared object file
#
# --no-build-isolation makes it compile against the pysam conda installed here
# instead. That requires the build dependencies (cython, numpy, setuptools, wheel,
# pysam) to already be present, which the environment yaml guarantees.
#
set -euo pipefail

pip install --no-build-isolation \
    "git+https://github.com/open2c/pairtools@fix-scaling"
