import glob, os, subprocess, sys

# this finds the packages to be excluded on the fly
# because windows is unpredictable
wheel, dest = sys.argv[1], sys.argv[2]
libs = r"C:\mm_env\Lib\site-packages\rdkit.libs"
names = [os.path.basename(p) for p in glob.glob(os.path.join(libs, "*.dll"))]
cmd = [
    "delvewheel", "repair",
    "--add-path", libs,
    "--no-dll", ";".join(names),
    "-w", dest, wheel,
]
print(" ".join(cmd))
subprocess.run(cmd, check=True)
