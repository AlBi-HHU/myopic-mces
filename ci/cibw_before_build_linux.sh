#!/bin/bash
# ai helped
set -euxo pipefail

# Setup python and rdkit paths
PYTAG=$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
ENV_PATH=/tmp/mm_env
rm -rf "$ENV_PATH"

RDKIT_SOURCE_DIR=/tmp/rdkit

# Setup conda env
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
export MAMBA_ROOT_PREFIX=/tmp/mamba_root

./bin/micromamba create -y -p "$ENV_PATH" -c conda-forge \
    "python=${PYTAG}" \
    "eigen" 

"$ENV_PATH/bin/python" -m pip install --upgrade pip
"$ENV_PATH/bin/python" -m pip install rdkit

# get RDKIT headers for bindings
cd /tmp
rm -rf "$RDKIT_SOURCE_DIR"
curl -L \
    -o rdkit.tar.gz \
    https://github.com/rdkit/rdkit/archive/refs/tags/Release_2026_03_5.tar.gz

mkdir "$RDKIT_SOURCE_DIR"
tar -xzf rdkit.tar.gz \
    --strip-components=1 \
    -C "$RDKIT_SOURCE_DIR"

cd "$RDKIT_SOURCE_DIR"

test -f "$RDKIT_SOURCE_DIR/Code/GraphMol/ROMol.h"
echo "RDKit headers:"
echo "  $RDKIT_SOURCE_DIR/Code"

# build boost from source
BOOST_VERSION=1.84.0
BOOST_ROOT=/tmp/boost_${BOOST_VERSION}
curl -L \
    https://archives.boost.io/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION//./_}.tar.gz \
    -o /tmp/boost.tar.gz

tar -xzf /tmp/boost.tar.gz -C /tmp
cd /tmp/boost_1_84_0
./bootstrap.sh \
    --prefix="$ENV_PATH" \
    --with-libraries=python,serialization,iostreams

./b2 \
    -j2 \
    variant=release \
    link=shared \
    runtime-link=shared \
    install

cmake -S /tmp/rdkit -B /tmp/rdkit-build \
    -DCMAKE_PREFIX_PATH=/tmp/mm_env \
    -DBoost_ROOT=/tmp/mm_env \
    -DBOOST_ROOT=/tmp/mm_env \
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DRDK_BUILD_CAIRO_SUPPORT=OFF \
    -DRDK_BUILD_COORDGEN_SUPPORT=OFF \
    -DRDK_BUILD_INCHI_SUPPORT=OFF \
    -DRDK_BUILD_MAEPARSER_SUPPORT=OFF \
    -DRDK_BUILD_AVALON_SUPPORT=OFF \
    -DRDK_BUILD_FREETYPE_SUPPORT=OFF \
    -DRDK_BUILD_PGSQL=OFF

export RDKit_INCLUDE_DIR="$RDKIT_SOURCE_DIR/Code"
export RDKit_LIBRARY_DIR="$ENV_PATH/lib"

cd /
"$ENV_PATH/bin/python" -c 'import rdkit; print(rdkit.__version__); print(rdkit.__file__)'

RDKIT_LIB_DIR="$ENV_PATH/lib/python${PYTAG}/site-packages/rdkit.libs"

ln -sf "$(echo "$RDKIT_LIB_DIR"/libRDKitGraphMol-*.so.1)" \
       "$ENV_PATH/lib/libRDKitGraphMol.so"

ln -sf "$(echo "$RDKIT_LIB_DIR"/libRDKitDataStructs-*.so.1)" \
       "$ENV_PATH/lib/libRDKitDataStructs.so"

ln -sf "$(echo "$RDKIT_LIB_DIR"/libRDKitSmilesParse-*.so.1)" \
       "$ENV_PATH/lib/libRDKitSmilesParse.so"
