# ExtXYZ.jl

![GitHub Workflow Status](https://img.shields.io/github/workflow/status/libAtoms/ExtXYZ.jl/CI) [![docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://libAtoms.github.io/ExtXYZ.jl/dev) [![docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://libatoms.github.io/ExtXYZ.jl/stable)

This package provides Julia bindings for the [extxyz](https://github.com/libAtoms/extxyz) C library which implements a parser and writer for the extended XYZ file format used in materials and molecular modelling, following the [specification](https://github.com/libAtoms/extxyz#extended-xyz-specification-and-parsing-tools) set out in the extxyz repo.

**Maintainer:** James Kermode ([@jameskermode](https://github.com/jameskermode)).

## Installation

This package is registered in the General registry, so installation of the latest stable release is as simple as pressing `]` to enter `pkg>` mode in the Julia REPL, and then entering:

```julia
pkg> add ExtXYZ
```

or for the development version:

```julia
pkg> dev https://github.com/libAtoms/ExtXYZ.jl
```

## Related packages

The [JuLIP.jl](https://github.com/JuliaMolSim/JuLIP.jl) package is an optional - but recommended - companion. JuLIP can use `ExtXYZ.jl` to read and write extended XYZ files to/from `JuLIP.Atoms` instances, using the functions `JuLIP.read_extxyz()` and `JuLIP.write_extxyz()`.

Please open issues/PRs here with suggestions of other packages it would be useful to provide interfaces to.

## Basic Usage

Four key functions are exported: `read_frame()` and `write_frame()` for reading and writing single configurations (snapshots), respectively, and `read_frames()` and `write_frames()` for reading and writing trajectories. All functions can work with string filenames, an open `Base.IO` instance or (intended primarily for internal use) a C `FILE*` pointer, stored as a `Ptr{Cvoid}` type.

```julia
using ExtXYZ

frame = read_frame("input.xyz")  # single atomic configuration, represented as a Dict{String}{Any}
write_frame("output.xyz", frame) # write a single frame to `output.xyz`. 

frame10 = read_frame("input.xyz", 10) # read a specific frame, counting from 1 for first frame in file

all_frames = read_frames("seq.xyz")  # read all frames, returns Vector{Dict{String}{Any}}
frames = read_frames("seq.xyz", 1:4) # specific range of frames

write_frames("output.xyz", frames, append=true) # append four frames to output
```

The function `iread_frames()` provides lazy file-reading using a `Channel`:

```julia
for frame in iread_frames("input.xyz")
    process(frame) # do something with each frame
do
```

`write_frames()` can also be used for asynchronous writing by passing in a `Channel`:

```julia
Channel() do ch
    @async write_frames(outfile, ch)
    
    for frame in frames
        put!(ch, frame)
    end
end
```

## Atoms data structure

In lieu of a package-independent data structure for representing atomic structures (i.e. an equivalent to ASE's `Atoms` class in the Python ecosystem), this package uses a `Dict{String}{Any}`. For the extended XYZ file:

```
8
Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.05.44" Properties=species:S:1:pos:R:3 Time=0.0
Si        0.00000000      0.00000000      0.00000000
Si        1.36000000      1.36000000      1.36000000
Si        2.72000000      2.72000000      0.00000000
Si        4.08000000      4.08000000      1.36000000
Si        2.72000000      0.00000000      2.72000000
Si        4.08000000      1.36000000      4.08000000
Si        0.00000000      2.72000000      2.72000000
Si        1.36000000      4.08000000      4.08000000
```

The internal representation, shown in JSON format for readability, is as follows:

```json
{
   "N_atoms": 8,
   "arrays": {
      "pos": [
         [
            0.0,
            0.0,
            0.0
         ],
         [
            1.36,
            1.36,
            1.36
         ],
         [
            2.72,
            2.72,
            0.0
         ],
         [
            4.08,
            4.08,
            1.36
         ],
         [
            2.72,
            0.0,
            2.72
         ],
         [
            4.08,
            1.36,
            4.08
         ],
         [
            0.0,
            2.72,
            2.72
         ],
         [
            1.36,
            4.08,
            4.08
         ]
      ],
      "species": [
         "Si",
         "Si",
         "Si",
         "Si",
         "Si",
         "Si",
         "Si",
         "Si"
      ]
   },
   "info": {
      "Time": 0.0
   },
   "cell": [
      [
         5.44,
         0.0,
         0.0
      ],
      [
         0.0,
         5.44,
         0.05
      ],
      [
         0.0,
         0.0,
         0.44
      ]
   ]
}
```

Important dictionary keys include:

 - `N_atoms` - the number of atoms (mandatory)
 - `cell` - the unit cell, a 3x3 matrix of floats containing the cell vectors as rows, i.e. the same as [ASE](https://wiki.fysik.dtu.dk/ase/ase/cell.html#ase.cell.Cell) (mandatory)
 - `pbc` - periodic boundary conditions, `Vector{Bool}` of length 3 (optional)
 - `info` - dictionary containing per-configuration key/value pairs parsed from the comment (line #2 in each frame). These can include scalars, vectors and matrices of integer, real, bool and string scalars or vectors. (mandatory, can be empty)
 - `arrays` - dictionary containing per-atom properties as a `N_component x N_atoms` matrix, reduced to a vector for the case `N_component = 1`. These represent scalar (`N_component = 1`) or vector (`N_component > 1`) per-atom properties, of integer (`I`), real (`R`), bool, (`L`) or string (`S`, scalars only) type. The set of properties is extracted from the special `Properties` key in the comment line. (mandatory, and must contain at least a string property `"species"` containing atomic symbols and a 3-column vector property


