# ExtXYZ.jl

This package provides Julia bindings for the [extxyz](https://github.com/libAtoms/extxyz) C library which implements a parser and writer for the extended XYZ file format used in materials and molecular modelling, according to the [specification](https://github.com/libAtoms/extxyz#extended-xyz-specification-and-parsing-tools) set out in the extxyz repo.

**Maintainer:** James Kermode ([@jameskermode](https://github.com/jameskermode)).

## Installation

This package will shortly be registered in the General registry, so installation of the latest stable release is as simple as:

```julia
] add ExtXYZ
```

or for the development version:

```julia
] dev https://github.com/libAtoms/ExtXYZ.jl
```

The [JuLIP.jl](https://github.com/JuliaMolSim/JuLIP.jl) package is an optional - but recommended - companion. JuLIP will shortly use `ExtXYZ.jl` internally to read and write extended XYZ files.

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

## Atoms data structure

In lieu of a package-independent data structure for representing atomic structures (i.e. an equivalent to ASE's `Atoms` class in the Python ecosystem), this package uses a `Dict{String}{Any}`. For an extended XYZ as follows

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

The interal representation, shown in JSON format for readability, is as follows:

```json
{
    "N_atoms": 8,
    "arrays": {
        "pos": [
            [
                0.0,
                0.0,
                1.36
            ],
            [
                0.0,
                4.08,
                4.08
            ],
            [
                0.0,
                4.08,
                0.0
            ],
            [
                1.36,
                1.36,
                2.72
            ],
            [
                1.36,
                2.72,
                2.72
            ],
            [
                1.36,
                0.0,
                1.36
            ],
            [
                2.72,
                2.72,
                4.08
            ],
            [
                2.72,
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
            0.0
        ],
        [
            0.0,
            0.0,
            5.44
        ]
    ]
}
```

Important top-level keys include:

 - `N_atoms` - the number of atoms (mandatory)
 - `cell` - the unit cell, a 3x3 matrix of floats containing the cell vectors as row, i.e. the same as [ASE](https://wiki.fysik.dtu.dk/ase/ase/cell.html#ase.cell.Cell) (mandatory)
 - `pbc` - periodic boundary conditions, `Vector{Bool}` of length 3 (optional)
 - `info` - dictionary containing per-configuration key/value pairs parsed from the comment (line #2 in each frame). These can include scalars, vectors and matrices of integer, real, bool and string scalars or vectors. (mandatory, can be empty)
 - `arrays` - dictionary containing per-atom properties as a `N_component x `N_atoms` matrix, reduced to a vector for the case `N_component = 1`. These represent scalar (`N_component = 1`) or vector (`N_component > 1`) per-atom properties, of integer (`I`), real (`R`), bool, (`L`) or string (`S`, scalars only) type. The set of properties is extracted from the special `Properties` key in the comment line. (mandatory, and must contain at least a string property `"species"` containing atomic symbols and a 3-column vector property

`JuLIP.XYZ.read_extxyz()` and `JuLIP.XYZ.write_extxyz()` contain functionality to convert these dictionaries to/from `JuLIP.Atoms` type instances.
