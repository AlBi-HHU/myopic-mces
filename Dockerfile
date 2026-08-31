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

COPY pyproject.toml uv.lock ./
RUN --mount=type=ssh uv sync --locked --extra db --extra hdf5 --no-install-project

COPY . .
RUN --mount=type=ssh uv sync --locked --extra db --extra hdf5

ENV PATH="/app/.venv/bin:$PATH"
