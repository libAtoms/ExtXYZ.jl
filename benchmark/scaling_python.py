"""Python extxyz scaling on the same files the Julia benchmark generated."""
import time
import extxyz
import ase.io
import ase_extxyz.io  # noqa: F401 - registers the cextxyz format

SIZES = [10, 100, 1000, 4000, 16000, 64000, 200000]
REPEATS = 3
DIR = "/tmp/extxyz_bench"


def best_of(fn, n=REPEATS):
    best = float("inf")
    for _ in range(n):
        t0 = time.perf_counter()
        r = fn()
        el = time.perf_counter() - t0
        if isinstance(r, list):
            _ = r[-1]
        best = min(best, el)
    return best


print("natoms,read_regex_s,read_tok_s,write_s,ase_read_s,ase_write_s")
for natoms in SIZES:
    f = f"{DIR}/bench_{natoms}.xyz"
    dicts = extxyz.read_dicts(f)
    t_read = best_of(lambda: extxyz.read_dicts(f))
    t_tok = best_of(lambda: extxyz.read_dicts(f, use_regex=False))
    t_write = best_of(lambda: extxyz.write_dicts(f"{DIR}/out_py_{natoms}.xyz", dicts))
    # ASE Atoms via the ase-extxyz plugin: Python analogue of AtomsBase load/save
    atoms = ase.io.read(f, format="cextxyz")
    t_ase_read = best_of(lambda: ase.io.read(f, format="cextxyz"))
    t_ase_write = best_of(
        lambda: ase.io.write(f"{DIR}/out_ase_{natoms}.xyz", atoms, format="cextxyz"))
    print(f"{natoms},{t_read:.6g},{t_tok:.6g},{t_write:.6g},{t_ase_read:.6g},{t_ase_write:.6g}")
