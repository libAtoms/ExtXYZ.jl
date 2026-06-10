# Manual benchmark for release notes: run once with the extxyz_jll 0.1.3
# artifact and once with libextxyz 0.4.0 (via ~/.julia/artifacts/Overrides.toml)
# and compare. Not part of the CI test suite.
#
#   julia --project=. benchmark/benchmark.jl [nframes] [natoms]

using ExtXYZ
using extxyz_jll
using Printf

nframes = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1000
natoms  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 64

println("libextxyz: ", extxyz_jll.libextxyz)
println("frames: $nframes, atoms/frame: $natoms")

frame = Dict{String,Any}(
    "N_atoms" => natoms,
    "cell" => [5.44 0.0 0.0; 0.0 5.44 0.0; 0.0 0.0 5.44],
    "pbc" => [true, true, true],
    "info" => Dict{String,Any}("energy" => -42.0, "config_type" => "bulk"),
    "arrays" => Dict{String,Any}(
        "species" => fill("Si", natoms),
        "pos" => rand(3, natoms),
        "forces" => rand(3, natoms)))
frames = [frame for _ in 1:nframes]

file = tempname() * ".xyz"
try
    # warm up, then time
    write_frames(file, frames[1:10])
    read_frames(file)

    t_write = @elapsed write_frames(file, frames)
    t_read  = @elapsed read_frames(file)
    t_tok   = try
        @elapsed read_frames(file; use_regex=false)
    catch
        NaN  # not supported by libextxyz < 0.4
    end

    @printf "write:            %8.3f s  (%8.1f frames/s)\n" t_write nframes/t_write
    @printf "read (regex):     %8.3f s  (%8.1f frames/s)\n" t_read nframes/t_read
    @printf "read (tokenizer): %8.3f s  (%8.1f frames/s)\n" t_tok nframes/t_tok
finally
    rm(file, force=true)
end
