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

# ---------------------------------------------------------------------------
# Stop pairtools leaking bgzip's version banner into its own stdout.
#
# pairtools/lib/fileio.py probes bgzip before using it:
#
#     subprocess.Popen(shlex.split('bgzip --version'),
#                      stderr=subprocess.PIPE, text=True)
#
# It captures stderr but not stdout, and bgzip prints
#
#     bgzip (htslib) 1.24
#     Copyright (C) 2026 Genome Research Ltd.
#
# to *stdout*, which the subprocess inherits from pairtools. So whenever a
# pairtools command reads a .gz input and writes pairs to stdout, those two
# lines land ahead of the pairs header and corrupt the stream. In this workflow
# that is `pairtools merge <run>.pairs.gz ... | pairtools dedup ...` in rule
# merge_dedup, which makes every library with more than one run fail with
# "Input file is not valid .pairs, has no header or is empty". Commands writing
# to a file with -o are unaffected -- there the banner only reaches the log,
# which is why single-run libraries work.
#
# This is a no-op if already patched or if the probe is gone, so it will not
# break once pairtools fixes this upstream.
# ---------------------------------------------------------------------------
python <<'PATCH_PY'
import pathlib

NEEDLE = ("subprocess.Popen(shlex.split('bgzip --version'), "
          "stderr=subprocess.PIPE, text=True)")
REPLACEMENT = ("subprocess.Popen(shlex.split('bgzip --version'), "
               "stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)")

try:
    import pairtools.lib.fileio as fileio
except ImportError as exc:
    print(f"bgzip-banner patch: pairtools not importable ({exc})")
else:
    path = pathlib.Path(fileio.__file__)
    source = path.read_text()
    if "stdout=subprocess.DEVNULL" in source:
        print("bgzip-banner patch: already applied")
    elif NEEDLE not in source:
        print(f"bgzip-banner patch: probe not found in {path}; upstream likely fixed")
    else:
        path.write_text(source.replace(NEEDLE, REPLACEMENT, 1))
        print(f"bgzip-banner patch: applied to {path}")
PATCH_PY
