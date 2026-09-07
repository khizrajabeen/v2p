# v2p — variant sequences to protein sequences.
#
# Two stages, on purpose. The references are ~18 GB unpacked and the genome
# alone is 3.1 GB; baking them into the default image would produce
# something nobody can pull, and would pin a reference release into every
# derived image. So:
#
#   docker build -t v2p .                            # ~150 MB, no references
#   docker build -t v2p:refs --target with-refs .    # large, references baked
#
# The reference *release* is pinned either way (the ARGs below), so a run is
# reproducible whichever image is used. The slim image expects the
# references mounted:
#
#   docker run --rm -v /data/refs:/refs -v $PWD/calls:/calls -v $PWD/out:/out \
#       v2p run /calls --ref /refs --outdir /out/release --logdir /out/logs
#
# Note --logdir outside --outdir: MANIFEST.txt is written as the last step
# of packaging, so anything written into the release afterwards makes it
# stale and the I8 invariant will correctly fail the run.

# --------------------------------------------------------------------------
FROM python:3.12-slim AS base

# The reference release, pinned. Changing these changes the output, so they
# are build arguments rather than something the entrypoint decides.
ARG GENCODE_RELEASE=44
ARG GENOME_BUILD=GRCh38
ENV V2P_GENCODE_RELEASE=${GENCODE_RELEASE} \
    V2P_GENOME_BUILD=${GENOME_BUILD} \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# pyfaidx needs no compiler. curl and make are for the reference fetch
# stage and the build assertion in scripts/00_fetch_references.sh.
RUN apt-get update \
 && apt-get install -y --no-install-recommends curl ca-certificates make \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/v2p
COPY pyproject.toml requirements.txt README.md ./
COPY src/ ./src/
COPY scripts/ ./scripts/
COPY config/ ./config/
COPY tests/ ./tests/
COPY docs/ ./docs/
COPY Makefile ./

RUN pip install --no-cache-dir .

# Prove the image works before it is tagged. This needs no network and no
# reference download, and it fails the build rather than shipping an image
# whose own tests do not pass.
RUN python tests/test_pipeline.py \
 && python tests/test_invariants.py \
 && python tests/test_config.py

# A pip-installed v2p cannot yet locate the numbered stage scripts, so the
# checkout stays on PATH and is what the entrypoint runs. See stage_dir()
# in src/v2p/cli.py for the search order and the error it raises.
ENV PATH="/opt/v2p/scripts:${PATH}"

WORKDIR /work
ENTRYPOINT ["v2p"]
CMD ["--help"]

# --------------------------------------------------------------------------
# The same image with the pinned references baked in. Large, and built only
# when asked for with --target with-refs.
FROM base AS with-refs

ARG GENCODE_RELEASE=44
RUN bash /opt/v2p/scripts/00_fetch_references.sh /refs \
 && rm -f /refs/*.fa.gz

# The audit checks that the GTF, genome and proteome are being *read*
# correctly, separately from whether the output looks right. Running it at
# build time means a reference release that changed its interface fails
# here rather than silently in someone's output.
RUN v2p audit --ref /refs

WORKDIR /work
ENTRYPOINT ["v2p"]
CMD ["--help"]
