# ai helped
$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

# Setup python and rdkit paths
$PYTAG        = & python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
$PYTAG_NO_DOT = & python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')"
$EnvPath = "C:\mm_env"
if (Test-Path $EnvPath) { Remove-Item -Recurse -Force $EnvPath }

$RdkitSourceDir = "C:\rdkit_src"

# Setup conda env
Invoke-WebRequest -Uri "https://micro.mamba.pm/api/micromamba/win-64/latest" -OutFile "C:\micromamba.tar.bz2"
# micromamba's win-64 archive layout: Library/bin/micromamba.exe
tar -xjf "C:\micromamba.tar.bz2" -C "C:\" "Library/bin/micromamba.exe"
$Micromamba = "C:\Library\bin\micromamba.exe"
$env:MAMBA_ROOT_PREFIX = "C:\mamba_root"

& $Micromamba create -y -p $EnvPath -c conda-forge "python=$PYTAG" "eigen"

& "$EnvPath\python.exe" -m pip install --upgrade pip
& "$EnvPath\python.exe" -m pip install "rdkit==2024.09.1"

$RDKitDir = & "$EnvPath\python.exe" -c "import rdkit; from pathlib import Path; print(Path(rdkit.__file__).parent)"
Write-Host $RDKitDir
Get-ChildItem "$RDKitDir\..\rdkit.libs"

# Compiler: rely on MSVC already provided by the cibuildwheel Windows runner image.
# Locate and load the VS dev environment (vcvarsall) so cl.exe / link.exe are on PATH.
$VsWhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$VsPath  = & $VsWhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
$VcVars  = Join-Path $VsPath "VC\Auxiliary\Build\vcvars64.bat"

# vcvarsall sets env vars in cmd's session, not PowerShell's — capture and re-apply
cmd.exe /c "call `"$VcVars`" && set" | ForEach-Object {
    if ($_ -match "^(.*?)=(.*)$") {
        [System.Environment]::SetEnvironmentVariable($matches[1], $matches[2])
    }
}
cl.exe 2>&1 | Select-String "Version"

# get RDKit headers for bindings
Set-Location "C:\"
if (Test-Path $RdkitSourceDir) { Remove-Item -Recurse -Force $RdkitSourceDir }
Invoke-WebRequest `
    -Uri "https://github.com/rdkit/rdkit/archive/refs/tags/Release_2024_09_1.tar.gz" `
    -OutFile "C:\rdkit.tar.gz"

New-Item -ItemType Directory -Path $RdkitSourceDir | Out-Null
tar -xzf "C:\rdkit.tar.gz" --strip-components=1 -C $RdkitSourceDir

Set-Location $RdkitSourceDir
if (-not (Test-Path "$RdkitSourceDir\Code\GraphMol\ROMol.h")) {
    throw "RDKit headers not found"
}
Write-Host "RDKit headers: $RdkitSourceDir\Code"

# determine Boost version used by pip RDKit
Set-Location "C:\"
$RDKitSite = & "$EnvPath\python.exe" -c "import rdkit, os; print(os.path.dirname(rdkit.__file__))"
$RDChem = Get-ChildItem -Path $RDKitSite -Filter "rdchem*.pyd" -Recurse | Select-Object -First 1

if (-not $RDChem) { throw "rdchem*.pyd not found" }

# NOTE: there's no direct Windows equivalent of readelf/otool for pulling a
# versioned soname the same way — DLLs on Windows don't embed a version in
# the filename the way libboost_python*.so.X.Y.Z does on Linux, and
# `dumpbin /dependents` (requires the VS dev shell, loaded above) just lists
# DLL names, not versions. In practice it's more reliable to pin the Boost
# version to match whatever conda-forge's rdkit==2024.09.1 build used, and
# verify against dumpbin output in CI logs rather than parsing it automatically.
$DependentsOutput = & dumpbin /dependents $RDChem.FullName
$BoostDllLine = $DependentsOutput | Select-String "boost_python" | Select-Object -First 1
Write-Host "Boost.Python DLL referenced: $BoostDllLine"

$BoostVersion = "1.85.0"   # verify/pin against actual RDKit conda-forge build
Write-Host "Using Boost version: $BoostVersion"

# build matching Boost from source
$BoostRoot = "C:\boost_$($BoostVersion -replace '\.','_')"
if (Test-Path $BoostRoot) { Remove-Item -Recurse -Force $BoostRoot }
$BoostUnderscored = $BoostVersion -replace '\.', '_'
Invoke-WebRequest `
    -Uri "https://archives.boost.io/release/$BoostVersion/source/boost_$BoostUnderscored.zip" `
    -OutFile "C:\boost.zip"
Expand-Archive -Path "C:\boost.zip" -DestinationPath "C:\boost_extract"
Move-Item "C:\boost_extract\boost_$BoostUnderscored" $BoostRoot

Set-Location $BoostRoot
& .\bootstrap.bat

& .\b2.exe `
    toolset=msvc `
    variant=release `
    link=shared `
    runtime-link=shared `
    --prefix="$EnvPath" `
    --with-python --with-serialization --with-iostreams --with-system `
    install

& "$EnvPath\python.exe" -c "import sys; print(sys.version)"

if (Test-Path "C:\rdkit-build") { Remove-Item -Recurse -Force "C:\rdkit-build" }

cmake -S "C:\rdkit_src" -B "C:\rdkit-build" `
    -DCMAKE_PREFIX_PATH="$EnvPath" `
    -DBoost_ROOT="$EnvPath" `
    -DBOOST_ROOT="$EnvPath" `
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF `
    -DRDK_BUILD_CPP_TESTS=OFF `
    -DRDK_BUILD_CAIRO_SUPPORT=OFF `
    -DRDK_BUILD_COORDGEN_SUPPORT=OFF `
    -DRDK_BUILD_INCHI_SUPPORT=OFF `
    -DRDK_BUILD_MAEPARSER_SUPPORT=OFF `
    -DRDK_BUILD_AVALON_SUPPORT=OFF `
    -DRDK_BUILD_FREETYPE_SUPPORT=OFF `
    -DRDK_BUILD_PGSQL=OFF

$env:RDKit_INCLUDE_DIR = "$RdkitSourceDir\Code"
$env:RDKit_LIBRARY_DIR = "$EnvPath\lib"

Set-Location "C:\"
& "$EnvPath\python.exe" -c "import rdkit; print(rdkit.__version__); print(rdkit.__file__)"

$env:RDKIT_LIB_DIR = "$EnvPath\Lib\site-packages\rdkit.libs"
$env:PATH = "$EnvPath\Lib\site-packages\rdkit.libs;$env:PATH"
$env:PYTAG = $PYTAG
$env:PYTAG_NO_DOT = $PYTAG_NO_DOT
