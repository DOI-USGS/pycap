import pathlib as pl
import subprocess

exampledir = pl.Path("../source/examples")
nbs = exampledir.glob("*.ipynb")
for nb in nbs:
    print("clearing", nb)
    cmd = (
        "jupyter",
        "nbconvert",
        "--ClearOutputPreprocessor.enabled=True",
        "--ClearMetadataPreprocessor.enabled=True",
        "--ClearMetadataPreprocessor.preserve_nb_metadata_mask={('kernelspec')}",
        "--inplace",
        nb,
    )
    proc = subprocess.run(cmd)
    assert proc.returncode == 0, f"Error running command: {' '.join(cmd)}"
    
    print("running fresh", nb)
    cmd = (
        "jupyter",
        "nbconvert",
        "--inplace",
        "--execute",
        nb,
    )
    proc = subprocess.run(cmd)
    assert proc.returncode == 0, f"Error running command: {' '.join(cmd)}"

