import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path):
    with open(path) as f:
        rows = [r for r in csv.DictReader(f)]
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


jl = load("/tmp/bench_julia_clean.csv")
py = load("/tmp/bench_python.csv")
n = jl["natoms"]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

ax = axes[0]
ax.loglog(n, jl["read_regex_s"], "o-", color="C2", label="Julia (regex)")
ax.loglog(n, jl["read_tok_s"], "o--", color="C2", alpha=0.6, label="Julia (tokenizer)")
ax.loglog(n, jl["load_s"], "d-", color="C4", label="Julia AtomsBase load (tokenizer)")
ax.loglog(n, py["read_regex_s"], "s-", color="C0", label="Python (regex)")
ax.loglog(n, py["read_tok_s"], "s--", color="C0", alpha=0.6, label="Python (tokenizer)")
ax.loglog(n, py["ase_read_s"], "v-", color="C5", label="Python ASE Atoms (tokenizer)")
ax.set_title("Read (1 frame)")
ax.set_xlabel("N atoms")
ax.set_ylabel("time (s)")
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)

ax = axes[1]
ax.loglog(n, jl["write_s"], "o-", color="C2", label="Julia")
ax.loglog(n, jl["save_s"], "d-", color="C4", label="Julia AtomsBase save")
ax.loglog(n, py["write_s"], "s-", color="C0", label="Python")
ax.loglog(n, py["ase_write_s"], "v-", color="C5", label="Python ASE Atoms")
ax.set_title("Write (1 frame)")
ax.set_xlabel("N atoms")
ax.set_ylabel("time (s)")
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)

ax = axes[2]
ax.loglog(n, jl["c_parse_s"], "^-", color="C3", label="C parse only")
ax.loglog(n, jl["dicts_s"], "^-", color="C1", label="+ dict conversion")
ax.loglog(n, jl["read_regex_s"], "o-", color="C2", label="full read_frames")
ax.loglog(n, jl["load_s"], "d-", color="C4", label="+ AtomsBase Atoms (load, tokenizer)")
ax.set_title("Julia read pipeline breakdown (regex)")
ax.set_xlabel("N atoms")
ax.set_ylabel("time (s)")
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)

fig.suptitle("ExtXYZ scaling: Julia (ExtXYZ.jl + libextxyz 0.4.0) vs Python (extxyz 0.4.1), same files, best of 3",
             fontsize=10)
fig.tight_layout()
out = "/Users/u1470235/.julia/dev/ExtXYZ/benchmark/scaling_julia_vs_python.png"
fig.savefig(out, dpi=150)
print(out)

# console summary at largest size
i = len(n) - 1
print(f"at N={int(n[i])}: read Jl/Py = {jl['read_regex_s'][i]/py['read_regex_s'][i]:.2f}, "
      f"tok Jl/Py = {jl['read_tok_s'][i]/py['read_tok_s'][i]:.2f}, "
      f"write Jl/Py = {jl['write_s'][i]/py['write_s'][i]:.2f}")
print(f"Julia overheads at N={int(n[i])}: convert = "
      f"{(jl['dicts_s'][i]-jl['c_parse_s'][i])/jl['dicts_s'][i]*100:.0f}% of dict read, "
      f"high-level+Channel = {(jl['read_regex_s'][i]-jl['dicts_s'][i])*1000:.1f} ms")
