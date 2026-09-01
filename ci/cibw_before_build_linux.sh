#!/bin/bash
# this is claude
set -euxo pipefail

PYTAG=$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
ENV_PATH=/tmp/mm_env
rm -rf "$ENV_PATH"

curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
export MAMBA_ROOT_PREFIX=/tmp/mamba_root

export CONDA_OVERRIDE_GLIBC="2.28"

./bin/micromamba create -y -p "$ENV_PATH" -c conda-forge \
    "python=${PYTAG}" \
    "boost=1.84" \
    "rdkit>=2022.09.5" \
    "librdkit>=2022.09.5" \
    "librdkit-dev>=2022.09.5" \
    "eigen" \
    "sysroot_linux-64=2.28"

# echo "--- Checking RDKit install ---"
# find /tmp/mm_env -iname "ROMol.h" -o -iname "libRDKit*" 2>/dev/null
# ls -la /tmp/mm_env/include/ 2>/dev/null || echo "include/ missing entirely"
# ls -la /tmp/mm_env/include/rdkit/ 2>/dev/null || echo "include/rdkit/ missing"
