"""
This package provides Julia bindings for the
[extxyz](https://github.com/libAtoms/extxyz) C library which implements a parser
and writer for the extended XYZ file format used in materials and molecular
modelling, following the
[specification](https://github.com/libAtoms/extxyz#extended-xyz-specification-and-parsing-tools)
set out in the extxyz repo.

## Basic Usage

Four key functions are exported: `read_frame()` and `write_frame()` for reading
and writing single configurations (snapshots), respectively, and `read_frames()`
and `write_frames()` for reading and writing trajectories. All functions can
work with string filenames, an open `Base.IO` instance or (intended primarily
for internal use) a C `FILE*` pointer, stored as a `Ptr{Cvoid}` type.

```julia
using ExtXYZ

frame = read_frame("input.xyz")  # single atomic configuration, represented as a Dict{String}{Any}
write_frame("output.xyz", frame) # write a single frame to `output.xyz`. 

frame10 = read_frame("input.xyz", 10) # read a specific frame, counting from 1 for first frame in file

all_frames = read_frames("seq.xyz")  # read all frames, returns Vector{Dict{String}{Any}}
frames = read_frames("seq.xyz", 1:4) # specific range of frames

write_frames("output.xyz", frames, append=true) # append four frames to output
```
"""
module ExtXYZ

include("fileio.jl") 
include("atoms.jl")

end
