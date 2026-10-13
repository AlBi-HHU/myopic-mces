import subprocess, sys, tempfile, venv, glob, os, sysconfig
from pathlib import Path

wheels_dir = Path(sys.argv[1]).resolve()
smoke_test = Path(sys.argv[2]).resolve()

# get wheel that fits python version
tag = f"cp{sys.version_info.major}{sys.version_info.minor}"
candidates = sorted(wheels_dir.glob(f"*-{tag}-{tag}-*.whl"))
if len(candidates) != 1:
    sys.exit(f"expecting exactly one wheel with {tag} in {wheels_dir}, found: {candidates}")
wheel = candidates[0]
print("--- Testing:", wheel.name)

with tempfile.TemporaryDirectory() as tmp:
    env = Path(tmp) / "venv"
    venv.EnvBuilder(with_pip=True, clear=True).create(env)
    py = env / ("Scripts/python.exe" if os.name == "nt" else "bin/python")

    # install exactly this wheel
    subprocess.run([str(py), "-m", "pip", "install", str(wheel)], check=True)

    # macOS: import the extension on its own first, it has to find RDKit's
    # libraries via its rpath without rdkit being imported beforehand
    # (on linux/windows rdkit has to be imported first, which the package does)
    if sys.platform == "darwin":
        subprocess.run([str(py), "-c", "import mmces_filters"], check=True, cwd=tmp)

    # run wheel in external tempfile to check
    subprocess.run([str(py), str(smoke_test)], check=True, cwd=tmp)
