#!/bin/bash
set -euo pipefail

PYTAG=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
ENV_PATH=/tmp/mm_env
rm -rf "$ENV_PATH"

curl -Ls https://micro.mamba.pm/api/micromamba/osx-$(uname -m | sed 's/x86_64/64/;s/arm64/arm64/')/latest | tar -xj bin/micromamba
export MAMBA_ROOT_PREFIX=/tmp/mamba_root

./bin/micromamba create -y -p "$ENV_PATH" -c conda-forge \
    "python=${PYTAG}" "boost=1.84" rdkit
