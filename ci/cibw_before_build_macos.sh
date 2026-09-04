#!/bin/bash
# ai helped
set -euxo pipefail

# Setup python and rdkit paths
PYTAG=$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PYTAG_NO_DOT=$(python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')")
ENV_PATH=/tmp/mm_env
rm -rf "$ENV_PATH"

RDKIT_SOURCE_DIR=/tmp/rdkit

# Detect arch -> micromamba platform tag
ARCH=$(uname -m)
if [ "$ARCH" = "arm64" ]; then
    MAMBA_PLATFORM="osx-arm64"
else
    MAMBA_PLATFORM="osx-64"
fi

# Setup conda env
curl -Ls "https://micro.mamba.pm/api/micromamba/${MAMBA_PLATFORM}/latest" | tar -xj bin/micromamba
export MAMBA_ROOT_PREFIX=/tmp/mamba_root

./bin/micromamba create -y -p "$ENV_PATH" -c conda-forge \
    "python=${PYTAG}" \
    "eigen"

"$ENV_PATH/bin/python" -m pip install --upgrade pip
"$ENV_PATH/bin/python" -m pip install "rdkit==2024.09.1"

RDKit_DIR=$("$ENV_PATH/bin/python" - <<'PY'
import rdkit
from pathlib import Path
print(Path(rdkit.__file__).parent)
PY
)

echo "$RDKit_DIR"
ls -lh "$RDKit_DIR/../rdkit.libs"

# Compiler: use Xcode command-line-tools clang (already on GH/GitLab macOS runners)
export CC=clang
export CXX=clang++
$CC --version
$CXX --version

# get RDKit headers for bindings
cd /tmp
rm -rf "$RDKIT_SOURCE_DIR"
curl -L \
    -o rdkit.tar.gz \
    https://github.com/rdkit/rdkit/archive/refs/tags/Release_2024_09_1.tar.gz

mkdir "$RDKIT_SOURCE_DIR"
tar -xzf rdkit.tar.gz \
    --strip-components=1 \
    -C "$RDKIT_SOURCE_DIR"

cd "$RDKIT_SOURCE_DIR"

test -f "$RDKIT_SOURCE_DIR/Code/GraphMol/ROMol.h"
echo "RDKit headers:"
echo "  $RDKIT_SOURCE_DIR/Code"

# determine Boost version used by pip RDKit
cd /

RDKIT_SITE=$("$ENV_PATH/bin/python" -c '
import rdkit
import os
print(os.path.dirname(rdkit.__file__))
')
RDCHEM=$(find "$RDKIT_SITE" -name 'rdchem*.so' -print -quit)
test -n "$RDCHEM"

# NOTE: macOS equivalent of readelf -d is otool -L (lists linked dylibs,
# not versioned .so-style sonames the same way Linux does). RDKit's macOS
# wheels typically link boost_python via an @rpath dylib name like
# libboost_python312.dylib or libboost_python-mt-x64.dylib — the exact
# naming has varied across RDKit/conda-forge releases, so verify this
# actually matches what you see in CI output and adjust the regex below.
BOOST_SONAME=$(otool -L "$RDCHEM" |
    sed -n 's/.*\/\(libboost_python[^ ]*\.dylib\).*/\1/p' |
    head -n1)

test -n "$BOOST_SONAME"
# Try to pull a version-looking token out of the name; if RDKit's dylib
# name doesn't embed a version (common on macOS/conda-forge), you'll need
# to instead pin BOOST_VERSION manually here based on what conda-forge's
# rdkit==2024.09.1 build actually depends on for boost.
BOOST_VERSION=$(echo "$BOOST_SONAME" |
    sed -n 's/.*-\([0-9]\+\.[0-9]\+\.[0-9]\+\)\.dylib$/\1/p')

if [ -z "$BOOST_VERSION" ]; then
    echo "WARNING: could not auto-detect Boost version from dylib name '$BOOST_SONAME'"
    echo "Falling back to a pinned version — verify this matches your RDKit build."
    BOOST_VERSION="1.85.0"
fi

echo "RDKit Python extension:"
echo "  $RDCHEM"
echo "RDKit Boost.Python library:"
echo "  $BOOST_SONAME"
echo "RDKit Boost version:"
echo "  $BOOST_VERSION"

# build matching Boost from source
BOOST_ROOT="/tmp/boost_${BOOST_VERSION}"
rm -rf "$BOOST_ROOT"
curl -L \
    "https://archives.boost.io/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION//./_}.tar.gz" \
    -o /tmp/boost.tar.gz
mkdir "$BOOST_ROOT"
tar -xzf /tmp/boost.tar.gz \
    --strip-components=1 \
    -C "$BOOST_ROOT"
cd "$BOOST_ROOT"
./bootstrap.sh \
    --prefix="$ENV_PATH" \
    --with-libraries=python,serialization,iostreams,system

./b2 \
    variant=release \
    link=shared \
    runtime-link=shared \
    install

"$ENV_PATH/bin/python" -c 'import sys; print(sys.version)'
otool -L "$ENV_PATH/lib/libboost_python${PYTAG_NO_DOT}.dylib" || true

rm -rf /tmp/rdkit-build

cmake -S /tmp/rdkit -B /tmp/rdkit-build \
    -DCMAKE_PREFIX_PATH="$ENV_PATH" \
    -DBoost_ROOT="$ENV_PATH" \
    -DBOOST_ROOT="$ENV_PATH" \
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

export RDKIT_LIB_DIR="$ENV_PATH/lib/python${PYTAG}/site-packages/rdkit.libs"
export DYLD_LIBRARY_PATH="$ENV_PATH/lib/python${PYTAG}/site-packages/rdkit.libs:${DYLD_LIBRARY_PATH:-}"
export PYTAG
export PYTAG_NO_DOT
