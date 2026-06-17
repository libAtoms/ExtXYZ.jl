# Scaling benchmark for ExtXYZ.jl mirroring upstream benchmarks/bench_read.py:
# sizes 10..200k atoms, 1 frame, Cu + forces + energy/step info, best-of-3.
# Also stages the read pipeline to isolate Julia-side overheads:
#   c_parse   - extxyz_read_ll_opts + free_dict only (no conversion)
#   dicts     - read_frame_dicts (C parse + C->Julia dict conversion)
#   full      - read_frames (adds high-level frame assembly + Channel)
using ExtXYZ
using extxyz_jll
using Printf

const SIZES = [10, 100, 1000, 4000, 16000, 64000, 200000]
const REPEATS = 3
const DIR = "/tmp/extxyz_bench"

function make_frame(natoms)
    a = 3.615  # Cu lattice constant
    n = ceil(Int, cbrt(natoms / 4))
    pos = Matrix{Float64}(undef, 3, natoms)
    k = 0
    basis = [(0.0,0.0,0.0), (0.5,0.5,0.0), (0.5,0.0,0.5), (0.0,0.5,0.5)]
    for i in 0:n-1, j in 0:n-1, l in 0:n-1, b in basis
        k += 1
        k > natoms && break
        pos[:, k] .= a .* (i + b[1], j + b[2], l + b[3])
        k == natoms && break
    end
    L = n * a
    Dict{String,Any}(
        "N_atoms" => natoms,
        "cell" => [L 0.0 0.0; 0.0 L 0.0; 0.0 0.0 L],
        "pbc" => [true, true, true],
        "info" => Dict{String,Any}("energy" => -1.234, "step" => 42),
        "arrays" => Dict{String,Any}(
            "species" => fill("Cu", natoms),
            "pos" => pos,
            "forces" => 0.05 .* randn(3, natoms)))
end

best_of(f, n=REPEATS) = minimum((f(); @elapsed f()) for _ in 1:n)

function c_parse_time(file; use_regex=true)
    # pure C parse: grammar+file -> C dicts -> free, no Julia conversion
    ExtXYZ.cfopen(file, "r") do fp
        nat = Ref{Cint}(0)
        info = Ref{Ptr{ExtXYZ.DictEntry}}()
        arrays = Ref{Ptr{ExtXYZ.DictEntry}}()
        err = zeros(UInt8, 1024)
        res = ccall((:extxyz_read_ll_opts, extxyz_jll.libextxyz), Cint,
                    (Ptr{Cvoid}, Ptr{Cvoid}, Ref{Cint}, Ptr{Ptr{ExtXYZ.DictEntry}},
                     Ptr{Ptr{ExtXYZ.DictEntry}}, Cstring, Ptr{UInt8}, Cint),
                    ExtXYZ._kv_grammar[], fp, nat, info, arrays, C_NULL, err,
                    use_regex ? 0 : 1)
        res == 1 || error("parse failed")
        ExtXYZ.cfree_dict(info[])
        ExtXYZ.cfree_dict(arrays[])
    end
end

mkpath(DIR)
println("natoms,file_mb,write_s,read_regex_s,read_tok_s,c_parse_s,dicts_s,load_s,save_s")
for natoms in SIZES
    frame = make_frame(natoms)
    file = joinpath(DIR, "bench_$natoms.xyz")
    t_write = best_of(() -> write_frames(file, [frame]))
    mb = filesize(file) / 1e6
    t_read  = best_of(() -> read_frames(file; use_regex=true))
    t_tok   = best_of(() -> read_frames(file; use_regex=false))
    t_c     = best_of(() -> c_parse_time(file))
    t_dicts = best_of(() -> ExtXYZ.cfopen(fp -> ExtXYZ.read_frame_dicts(fp; use_regex=true), file, "r"))
    t_load  = best_of(() -> ExtXYZ.load(file))
    sys     = ExtXYZ.load(file, 1)
    outfile = joinpath(DIR, "out_atoms_$natoms.xyz")
    t_save  = best_of(() -> ExtXYZ.save(outfile, sys))
    @printf "%d,%.6f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n" natoms mb t_write t_read t_tok t_c t_dicts t_load t_save
end
