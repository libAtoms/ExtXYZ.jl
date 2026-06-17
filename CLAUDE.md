# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

ExtXYZ.jl provides Julia bindings for the [extxyz](https://github.com/libAtoms/extxyz) C library (shipped as `extxyz_jll`), which parses and writes the extended XYZ file format used in materials/molecular modelling. It also exposes parsed configurations as `ExtXYZ.Atoms`, an implementation of the [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) (v0.5) `AbstractSystem` interface.

## Commands

```bash
# Run the full test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# Run a single test file (tests have their own deps in test/Project.toml)
julia --project=. -e 'using Pkg, TestEnv; TestEnv.activate(); include("test/dict.jl")'

# Build docs locally / run doctests
julia --project=docs docs/make.jl
```

There are only two test files: `test/dict.jl` (Dict-based read/write round-trips) and `test/atomsbase.jl` (AtomsBase interface, uses `AtomsBaseTesting`). Both are included from `test/runtests.jl`.

CI (`.github/workflows/CI.yml`) tests Julia 1.10 and nightly on Linux/macOS/Windows, and runs doctests via Documenter.

## Architecture

Two source files, two layers:

**`src/fileio.jl` — C interop and file I/O.** Everything talks to `libextxyz` (>= 0.4.0) through `ccall`. Key pieces:
- `cfopen` overloads convert a filename, `IOStream`, or `IOBuffer` into a C `FILE*` (`Ptr{Cvoid}`), which is what the C library consumes. On Windows, filename open/close routes through `extxyz_fopen`/`extxyz_fclose` so the `FILE*` stays within the library's C runtime. All public read/write functions accept any of these. `IOBuffer` reading uses `fmemopen` and is **not supported on Windows**.
- `DictEntry` is a Julia mirror of the C library's linked-list dict struct. `convert` methods translate in both directions between `Ptr{DictEntry}` and the Julia `Dict{String,Any}` frame representation. ABI subtlety: per-atom string columns read from the C library arrive as one contiguous fixed-width NUL-padded buffer (`n_in_row < 0`, cell width `-n_in_row`), while Julia-constructed dicts use the legacy `char**` layout marked by `n_in_row = 0` — both the reader and the C `free_dict` branch on the sign.
- Reads go through `extxyz_read_ll_opts`, which requires a real (non-NULL) 1024-byte error buffer; `_is_eof_message` distinguishes end-of-file (empty message or all-whitespace natoms line) from genuine parse errors, and `iread_frames` relies on `EOFError` for termination.
- `__init__` compiles the extxyz key/value grammar once (`_kv_grammar`) for the C parser.
- Exports: `read_frame`, `read_frames`, `iread_frames` (lazy, Channel-based), `write_frame`, `write_frames`. `write_frames` also accepts a `Channel` for asynchronous writing. Reads take `use_regex::Bool=false` (default = fast whitespace tokenizer; `true` = stricter PCRE2 regex parser); writes take `fmt_i`/`fmt_f`/`fmt_b`/`fmt_s` printf-style format overrides via `extxyz_write_ll_fmt`.

**`src/atoms.jl` — AtomsBase layer.** Defines `Atoms{P,Q} <: AbstractSystem{3}` holding two NamedTuples (`atom_data`, `system_data`). Constructors convert from any `AbstractSystem` or from the frame `Dict`; `write_dict` goes back. `ExtXYZ.load`/`ExtXYZ.save` combine the two layers to read/write `Atoms` directly from/to files. Units follow ASE conventions (Å, eV, u; velocity unit is `sqrt(u"eV"/u"u")`). Only types the C layer can represent (integers, floats, strings, and arrays thereof) survive a round-trip; other properties are dropped with a warning. Spatial dimension is hard-coded to 3.

**Frame `Dict` representation** (the package-level data structure, see README for full detail):
- `"N_atoms"` — atom count
- `"cell"` — 3×3 matrix, cell vectors as rows (ASE convention)
- `"pbc"` — optional `Vector{Bool}` of length 3
- `"info"` — per-configuration key/value pairs from the comment line
- `"arrays"` — per-atom properties as `N_component × N_atoms` matrices (vectors when `N_component == 1`); must contain `"species"` and `"pos"`

The cell is stored under the `"Lattice"` key in the raw extxyz comment line; `extract_lattice!` moves it from `info` to the top-level `"cell"` key on read.
