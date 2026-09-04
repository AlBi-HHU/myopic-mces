#!/bin/bash
# ai helped
set -euxo pipefail

# Setup python and rdkit paths
PYTAG=$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PYTAG_NO_DOT=$(python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')")
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
"$ENV_PATH/bin/python" -m pip install "rdkit==2024.09.1"

RDKit_DIR=$("$ENV_PATH/bin/python" - <<'PY'
import rdkit
from pathlib import Path
print(Path(rdkit.__file__).parent)
PY
)

echo "$RDKit_DIR"
ls -lh "$RDKit_DIR/../rdkit.libs"


# install gcc 10 c++ compiler because otherwise we get glibc mismatch even though we are on the correct manylinux image??
dnf install -y gcc-toolset-10-gcc-c++
export CC=/opt/rh/gcc-toolset-10/root/usr/bin/gcc
export CXX=/opt/rh/gcc-toolset-10/root/usr/bin/g++
$CC --version
$CXX --version

# get RDKIT headers for bindings
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
BOOST_SONAME=$(readelf -d "$RDCHEM" |
    sed -n 's/.*Shared library: \[\(libboost_python[^]]*\)\].*/\1/p' |
    head -n1)

test -n "$BOOST_SONAME"
BOOST_VERSION=$(echo "$BOOST_SONAME" |
    sed -n 's/.*\.so\.\([0-9]\+\.[0-9]\+\.[0-9]\+\)$/\1/p')
test -n "$BOOST_VERSION"
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
    https://archives.boost.io/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION//./_}.tar.gz \
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
    -j2 \
    variant=release \
    link=shared \
    runtime-link=shared \
    install

"$ENV_PATH/bin/python" -c 'import sys; print(sys.version)'
strings "$ENV_PATH/lib/libboost_python${PYTAG_NO_DOT}.so" | grep GLIBCXX_ | sort -V | tail

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
export LD_LIBRARY_PATH=/tmp/mm_env/lib/python${PYTAG}/site-packages/rdkit.libs:$LD_LIBRARY_PATH
export PYTAG
export PYTAG_NO_DOT
