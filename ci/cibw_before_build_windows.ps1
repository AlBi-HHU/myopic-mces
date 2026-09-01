$ErrorActionPreference = "Stop"
$PyTag = (python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
$EnvPath = "C:\mm_env"
Remove-Item -Recurse -Force $EnvPath -ErrorAction SilentlyContinue

Invoke-WebRequest -Uri "https://micro.mamba.pm/api/micromamba/win-64/latest" -OutFile micromamba.tar.bz2
tar -xf micromamba.tar.bz2 Library/bin/micromamba.exe
$env:MAMBA_ROOT_PREFIX = "C:\mamba_root"

.\Library\bin\micromamba.exe create -y -p $EnvPath -c conda-forge python=$PyTag boost=1.84 rdkit
