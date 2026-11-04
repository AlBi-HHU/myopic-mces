# syntax=docker/dockerfile:1
FROM python:3.12-slim
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# git+ssh access to the private mcesdb repo at build time (forwarded via --ssh,
# never baked into the image). Build with:
#   DOCKER_BUILDKIT=1 docker build --ssh default -t myopic-mces .
# (needs the mces-database deploy key loaded in your local ssh-agent)
RUN apt-get update && apt-get install -y --no-install-recommends git openssh-client \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p -m 0700 /root/.ssh \
    && ssh-keyscan git.uni-jena.de >> /root/.ssh/known_hosts

WORKDIR /app
ENV UV_LINK_MODE=copy

COPY . .
RUN --mount=type=ssh --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --extra db --extra hdf5 --extra cplex

# `docplex config --upgrade` (see nextflow/modules.nf) switches the pip
# `cplex` package from its bundled size-limited Community Edition engine to a
# locally licensed CPLEX install at runtime, but only if the pip package's
# version matches the installed CPLEX Optimization Studio's version exactly.
# pyproject.toml/uv.lock only pin a floor for portability, so pin the exact
# version here to match this deployment's CPLEX Studio (bump when it's
# upgraded).
ARG CPLEX_PY_VERSION=22.2.0.1
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install --python /app/.venv/bin/python "cplex==${CPLEX_PY_VERSION}"

# Singularity/Apptainer runs task processes as the host user, not root, but
# this venv is built here as root. `docplex config --upgrade` writes into
# several spots under it (site-packages/cplex, and executables/libs under
# bin/), so make the whole venv writable rather than chase each one.
RUN chmod -R o+rwX /app/.venv

ENV PATH="/app/.venv/bin:$PATH"
