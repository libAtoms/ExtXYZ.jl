# ExtXYZ.jl

This package provides Julia bindings for the [extxyz](https://github.com/libAtoms/extxyz) C library which implements a parser and writer for the extended XYZ file format used in materials and molecular modelling, according to the [specification](https://github.com/libAtoms/extxyz#extended-xyz-specification-and-parsing-tools) set out in the extxyz repo.

**Maintainer:** James Kermode ([@jameskermode](https://github.com/jameskermode)).

## Installation

This package will shortly be registered in the General registry, so installation is as simple as:

```julia
] add ExtXYZ
```

for the development version:

```julia
] dev https://github.com/libAtoms/ExtXYZ.jl
```

The [JuLIP.jl](https://github.com/JuliaMolSim/JuLIP.jl) package is an optional - but recommended - companion. JuLIP will shortly use `ExtXYZ.jl` internally

## Basic Usage

Four main functions are exported: `read_frame()` and `write_frame()` for reading and writing single configurations (snapshots), respectively, and `read_frames()` and `write_frames()` for reading and writing trajectories.

```julia
using ExtXYZ

frame = read_frame("input.xyz")  # single atomic configuration, represented as a Dict{String}{Any}
write_frame("output.xyz", frame) # write a single frame to `output.xyz`

frame10 = read_frame("input.xyz", 10) # read a specific frame, counting from 1 for first frame in file

all_frames = read_frames("seq.xyz")  # read all frames
frames = read_frames("seq.xyz", 1:4) # specific range of frames
```

