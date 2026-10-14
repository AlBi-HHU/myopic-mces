# ai helped
$ErrorActionPreference = "Stop"
$ProgressPreference = 'SilentlyContinue'
Set-StrictMode -Version Latest

# Setup python and rdkit paths
$PYTAG        = & python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
$PYTAG_NO_DOT = & python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')"
$EnvPath = "C:\mm_env"
# unlink the Boost header junction first (see below), Remove-Item -Recurse
# would otherwise follow it and delete the cached Boost install
Get-ChildItem "$EnvPath\include" -Directory -Filter "boost-*" -ErrorAction SilentlyContinue |
    ForEach-Object { cmd.exe /c rmdir "$($_.FullName)" }
if (Test-Path $EnvPath) { Remove-Item -Recurse -Force $EnvPath }

$RdkitSourceDir = "C:\rdkit-src"

#setup conda env
$MicromambaDir = "C:\micromamba"

New-Item -ItemType Directory -Path $MicromambaDir -Force | Out-Null

Invoke-WebRequest `
    -Uri "https://micro.mamba.pm/api/micromamba/win-64/latest" `
    -OutFile "C:\micromamba.tar.bz2"

tar -xjf "C:\micromamba.tar.bz2" `
    -C $MicromambaDir `
    "Library/bin/micromamba.exe"

$Micromamba = "$MicromambaDir\Library\bin\micromamba.exe"
$env:MAMBA_ROOT_PREFIX = "C:\mamba_root"

& $Micromamba create -y -p $EnvPath -c conda-forge "python=$PYTAG" "eigen"

& "$EnvPath\python.exe" -m pip install --upgrade pip
& "$EnvPath\python.exe" -m pip install "rdkit==2025.3.4"

$RDKitDir = & "$EnvPath\python.exe" -c "import rdkit; from pathlib import Path; print(Path(rdkit.__file__).parent)"
Write-Host $RDKitDir
Get-ChildItem "$RDKitDir\..\rdkit.libs"

####### setup MSVC env
Write-Host "------------ setting up MSVC env"

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

Write-Host "------------ geting rdkit headers"
# get RDKit headers for bindings
Set-Location "C:\"
if (Test-Path $RdkitSourceDir) { Remove-Item -Recurse -Force $RdkitSourceDir }
Invoke-WebRequest `
    -Uri "https://github.com/rdkit/rdkit/archive/refs/tags/Release_2025_03_4.tar.gz" `
    -OutFile "C:\rdkit.tar.gz"

New-Item -ItemType Directory -Path $RdkitSourceDir | Out-Null
tar -xzf "C:\rdkit.tar.gz" --strip-components=1 -C $RdkitSourceDir

Set-Location $RdkitSourceDir
if (-not (Test-Path "$RdkitSourceDir\Code\GraphMol\ROMol.h")) {
    throw "RDKit headers not found"
}
Write-Host "RDKit headers: $RdkitSourceDir\Code"

$libsDir = "$EnvPath\Lib\site-packages\rdkit.libs"
foreach ($pat in "RDKitGraphMol*","RDKitDataStructs*","RDKitSmilesParse*","RDKitRDGeneral*","boost_python*") {
    $dll = Get-ChildItem $libsDir -Filter "$pat.dll" | Select-Object -First 1
    if (-not $dll) { throw "DLL $pat not found" }
    $exports = & dumpbin /exports $dll.FullName |
        Select-String '^\s+\d+\s+[0-9A-F]+\s+[0-9A-F]+\s+(\S+)' |
        ForEach-Object { $_.Matches[0].Groups[1].Value }
    $def = Join-Path $libsDir "$($dll.BaseName).def"
    "LIBRARY $($dll.Name)`r`nEXPORTS`r`n" + ($exports -join "`r`n") | Set-Content $def
    & lib /nologo /def:$def /out:"$libsDir\$($dll.BaseName).lib" /machine:x64
    if (-not (Test-Path "$libsDir\$($dll.BaseName).lib")) { throw "import lib for $($dll.Name) failed" }
}

# determine Boost version used by pip RDKit
Write-Host "----------- determining boost version"
Set-Location "C:\"
$RDKitSite = & "$EnvPath\python.exe" -c "import rdkit, os; print(os.path.dirname(rdkit.__file__))"
$RDChem = Get-ChildItem -Path $RDKitSite -Filter "rdchem*.pyd" -Recurse | Select-Object -First 1

if (-not $RDChem) { throw "rdchem*.pyd not found" }

# NOTE: there's no direct Windows equivalent of readelf/otool for pulling a
# versioned soname the same way — DLLs on Windows don't embed a version in
# the filename the way libboost_python*.so.X.Y.Z does on Linux, and
# `dumpbin /dependents` (requires the VS dev shell, loaded above) just lists
# DLL names, not versions. In practice it's more reliable to pin the Boost
# version to match whatever conda-forge's rdkit==2025.3.4 build used, and
# verify against dumpbin output in CI logs rather than parsing it automatically.
$DependentsOutput = & dumpbin /dependents $RDChem.FullName
$BoostDllLine = $DependentsOutput | Select-String "boost_python" | Select-Object -First 1
Write-Host "Boost.Python DLL referenced: $BoostDllLine"

$BoostVersion = "1.85.0"   # verify/pin against actual RDKit conda-forge build
Write-Host "Using Boost version: $BoostVersion"

# build matching Boost from source
# Boost is independent of the Python version (only headers + cmake configs are
# needed), so install it once into its own prefix and reuse it for all wheels
Write-Host "----------- building boost from source"
$BoostPrefix = "C:\boost-install"
$BoostUnderscored = $BoostVersion -replace '\.', '_'
$BoostRoot = "C:\boost_$BoostUnderscored"

if (Test-Path "$BoostPrefix\include\boost-*\boost\python.hpp") {
    Write-Host "Boost already installed in $BoostPrefix, skipping"
} else {
    if (Test-Path $BoostRoot) { Remove-Item -Recurse -Force $BoostRoot }

    # Expand-Archive on the zip takes ~30 min (tens of thousands of files),
    # the bundled bsdtar extracts the tarball in a fraction of that
    curl.exe -sSfL -o "C:\boost.tar.gz" `
        "https://archives.boost.io/release/$BoostVersion/source/boost_$BoostUnderscored.tar.gz"
    if ($LASTEXITCODE -ne 0) { throw "Boost download failed" }
    tar -xzf "C:\boost.tar.gz" -C "C:\"
    if ($LASTEXITCODE -ne 0) { throw "Boost extraction failed" }

    Set-Location $BoostRoot

    Write-Host "cl.exe:"
    where.exe cl.exe
    Write-Host "VSINSTALLDIR: $env:VSINSTALLDIR"
    Write-Host "VCToolsInstallDir: $env:VCToolsInstallDir"

    $ShimDir = "C:\boost-shim"
    New-Item -ItemType Directory -Force $ShimDir | Out-Null
    "@echo off`r`ncl.exe %*" | Set-Content "$ShimDir\msvc.bat"
    $env:Path = "$ShimDir;$env:Path"

    Write-Host "----------- bootstrapping Boost"
    & .\bootstrap.bat vc143 `
        --prefix="$BoostPrefix" `
        --with-libraries=python,serialization,iostreams,system
    if (-not (Test-Path ".\b2.exe")) { throw "Boost bootstrap failed" }

    # -d0 silences the ~17k "common.copy" lines of the header install
    & .\b2.exe `
        -j $env:NUMBER_OF_PROCESSORS -d0 `
        toolset=msvc address-model=64 variant=release link=shared runtime-link=shared `
        --prefix="$BoostPrefix" `
        --with-serialization --with-iostreams --with-system `
        install

    $toolsetFile = "$BoostPrefix\lib\cmake\BoostDetectToolset-1.85.0.cmake"
    (Get-Content $toolsetFile) -replace `
        '\(MSVC_VERSION GREATER 1929\) AND \(MSVC_VERSION LESS 1940\)', `
        '(MSVC_VERSION GREATER 1929) AND (MSVC_VERSION LESS 1960)' |
        Set-Content $toolsetFile
    Write-Host "=== BoostDetectToolset content (after patch) ==="
    Get-Content $toolsetFile | Select-String "vc143"

    Set-Location "C:\"
}

# the wheel build expects the headers in the conda env (BOOST_INCLUDE_DIR in
# pyproject.toml), which is recreated per Python version -> link them in
$BoostCachedInclude = (Get-ChildItem "$BoostPrefix\include" -Directory -Filter "boost-*" | Select-Object -First 1)
if (-not $BoostCachedInclude) { throw "Boost versioned include dir not found" }
New-Item -ItemType Directory -Force "$EnvPath\include" | Out-Null
$BoostIncludeDir = "$EnvPath\include\$($BoostCachedInclude.Name)"
New-Item -ItemType Junction -Path $BoostIncludeDir -Target $BoostCachedInclude.FullName | Out-Null
Write-Host "Boost include directory: $BoostIncludeDir"
$env:BOOST_INCLUDE_DIR = $BoostIncludeDir

Write-Host "----------- building rdkit"

if (Test-Path "C:\rdkit-build") {
    Remove-Item -Recurse -Force "C:\rdkit-build"
}

#-G "Ninja" `
cmake -S "C:\rdkit-src" -B "C:\rdkit-build" `
    -DCMAKE_PREFIX_PATH="$BoostPrefix;$EnvPath" `
    -DBoost_ROOT="$BoostPrefix" `
    -DBOOST_ROOT="$BoostPrefix" `
    -DEIGEN3_INCLUDE_DIR="$EnvPath\Library\include\eigen3" `
    -DEIGEN3_VERSION_OK=TRUE `
    -DBUILD_SHARED_LIBS=ON `
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF `
    -DRDK_BUILD_CPP_TESTS=OFF `
    -DRDK_BUILD_CAIRO_SUPPORT=OFF `
    -DRDK_BUILD_COORDGEN_SUPPORT=OFF `
    -DRDK_BUILD_INCHI_SUPPORT=OFF `
    -DRDK_USE_BOOST_SERIALIZATION=OFF `
    -DRDK_USE_BOOST_IOSTREAMS=OFF `
    -DRDK_BUILD_MAEPARSER_SUPPORT=OFF `
    -DRDK_BUILD_AVALON_SUPPORT=OFF `
    -DRDK_BUILD_FREETYPE_SUPPORT=OFF `
    -DRDK_BUILD_PGSQL=OFF

cmake --build "C:\rdkit-build" --target RDGeneral --config Release

Write-Host "=== RDKit generated headers ==="
Get-ChildItem "C:\rdkit-build" -Filter "export.h" -Recurse -ErrorAction SilentlyContinue | Select-Object FullName
Get-ChildItem "C:\rdkit-src" -Filter "export.h" -Recurse -ErrorAction SilentlyContinue | Select-Object FullName

Write-Host "=== Boost headers ==="
Test-Path "$BoostIncludeDir\boost*\python.hpp"

Write-Host "=== Boost libraries ==="
Get-ChildItem "$BoostPrefix\lib" -Recurse | Where-Object Name -match "boost" | Select-Object FullName

Set-Location "C:\"
& "$EnvPath\python.exe" -c "import rdkit; print(rdkit.__version__); print(rdkit.__file__)"

$env:RDKIT_LIB_DIR = "$EnvPath\Lib\site-packages\rdkit.libs"
$env:PATH = "$EnvPath\Lib\site-packages\rdkit.libs;$env:PATH"
$env:PYTAG = $PYTAG
$env:PYTAG_NO_DOT = $PYTAG_NO_DOT

Write-Host "=== RDKit ==="
Test-Path "C:\rdkit-src\Code\GraphMol\ROMol.h"
Test-Path "C:\rdkit-src\Code\RDGeneral\export.h"
Test-Path "C:\rdkit-build\Code\RDGeneral\export.h"

Write-Host "=== Boost ==="
Test-Path "$BoostIncludeDir\boost*\python.hpp"
Get-ChildItem "$BoostIncludeDir" -Filter "*boost*python*" -ErrorAction SilentlyContinue

if (-not (Test-Path "$BoostIncludeDir\boost\python.hpp")) {
    throw "Boost.Python headers were not installed"
}
