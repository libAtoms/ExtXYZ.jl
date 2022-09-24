using AtomsBase
using PeriodicTable
using StaticArrays
using Unitful

import Base: keys, getproperty, getindex, getfield, show, length, size, isapprox
import AtomsBase: bounding_box, boundary_conditions, species_type, position, atomic_symbol, atomic_number, atomic_mass, velocity

const D = 3 # TODO generalise to arbitrary spatial dimensions
const REQUIRED_PER_ATOM_SYMS = [:positions, :atomic_symbols, :atomic_numbers]
const REQUIRED_PER_SYSTEM_SYMS = [:box, :boundary_conditions]

const NAME_MAP = Base.ImmutableDict("pos" => "positions",
                                    "species" => "atomic_symbols",
                                     "Z" => "atomic_numbers",
                                     "mass" => "atomic_masses")
const REV_NAME_MAP = Base.ImmutableDict([value => key for (key, value) in NAME_MAP]...)

export Atoms

"""
`struct Atoms` is the main type for flexible systems of atoms
"""
# TODO: add paramteric entries for length and mass units
struct Atoms{P <: NamedTuple, Q <: NamedTuple} <: AbstractSystem{D} 
    atom_data::P
    system_data::Q

    # standard constructor - SYMS are the same 
    function Atoms{P, Q}(p::NT1, q::NT2) where {P <: NamedTuple{PSYMS}, NT1 <: NamedTuple{PSYMS}, 
                                                Q <: NamedTuple{QSYMS}, NT2 <: NamedTuple{QSYMS}} where {PSYMS, QSYMS}
        for sym in REQUIRED_PER_ATOM_SYMS
            sym ∉ PSYMS && error("Required per-atom symbol '$sym' missing in call to constructor")
            # TODO: add type checking for required symbols?
        end
        for sym in REQUIRED_PER_SYSTEM_SYMS
            sym ∉ QSYMS && error("Required per-system symbol '$sym' missing in call to constructor")
            # TODO: add type checking for required symbols?
        end
        new{NT1,NT2}(p, q)
    end
end

# outward facing constructor
Atoms(atom_data::NT1, system_data::NT2) where {NT1 <: NamedTuple, NT2 <: NamedTuple} = Atoms{NT1, NT2}(atom_data, system_data)

# ----------- conversion to/from dictionaries

_read_convert(value) = value
_read_convert(value::Int32) = Int(value)
_read_convert(value::Vector{Int32}) = Int.(value)

_write_convert(value) = value
_write_convert(value::Vector{T}) where {T<:Unitful.Length} = ustrip.(uconvert.(u"Å", value)) 
_write_convert(value::Vector{Vector{T}}) where {T<:Unitful.Length} = ustrip.(uconvert.(u"Å",(hcat((value)...))))
_write_convert(value::Vector{T}) where {T<:Unitful.Velocity} = ustrip.(uconvert.(u"Å/(Å*sqrt(u/eV))", value))
_write_convert(value::Vector{Vector{T}}) where {T<:Unitful.Velocity} = ustrip.(uconvert.(u"Å/(Å*sqrt(u/eV))",(hcat((value)...))))
_write_convert(value::Vector{T}) where {T<:Unitful.Mass} = ustrip.(uconvert.(u"u", value))
_write_convert(value::Vector{Vector{T}}) where {T<:Unitful.Mass} = ustrip.(uconvert.(u"u",(hcat((value)...))))
_write_convert(value::Vector{T}) where {T<:Unitful.Energy} = ustrip.(uconvert.(u"eV", value))
_write_convert(value::Vector{Vector{T}}) where {T<:Unitful.Energy} = ustrip.(uconvert.(u"eV",(hcat((value)...))))
_write_convert(value::Vector{T}) where {T} = ustrip.(value)
_write_convert(value::Vector{Vector{T}}) where {T} = ustrip.(hcat((value)...))
_write_convert(value::Vector{Symbol}) = string.(value)
_write_convert(value::Symbol) = string(value)

_dict_remap_names(d, names, pre, post, conv) = Dict{post}{Any}(post(get(names, pre(key), pre(key))) => conv(value) for (key, value) in pairs(d))
_dict_remap_fwd(d) = _dict_remap_names(d, NAME_MAP, identity, Symbol, _read_convert)
_dict_remap_rev(d) = _dict_remap_names(d, REV_NAME_MAP, string, String, _write_convert)

# read from dictionary of results as returned by read_frame()
function Atoms(dict::Dict{String}{Any})
    atom_data = _dict_remap_fwd(dict["arrays"])
    natoms = dict["N_atoms"]
    atom_data[:positions] = [ atom_data[:positions][:, i] for i=1:natoms ].*u"Å" # from matrix to vector of vectors

    if :atomic_symbols in keys(atom_data)
       sym = atom_data[:atomic_symbols] = Symbol.(atom_data[:atomic_symbols])
       Z = getfield.(elements[sym], :number)
       if haskey(atom_data, :atomic_numbers)
          all(atom_data[:atomic_numbers] .== Z) || error("inconsistent 'Z' and 'species' properties")
       else
          atom_data[:atomic_numbers] = Z
       end
    end
    atom_data[:atomic_numbers] === nothing && error("atomic numbers not defined - either 'Z' or 'species' must be present")
 
    # mass - lookup from atomic number if not present look up from atomic_symbols
    if !haskey(atom_data, :atomic_masses)
        if haskey(atom_data, :atomic_numbers)
            atom_data[:atomic_masses] = getfield.(elements[atom_data[:atomic_numbers]], :atomic_mass)
        elseif haskey(atom_data, :atomic_symbols)
            atom_data[:atomic_masses] = getfield.(elements[atom_data[:atomic_symbols]], :atomic_mass)
        end
    end
  
    system_data = _dict_remap_fwd(dict["info"])
    system_data[:box] = [dict["cell"][i, :] for i in 1:D ].*u"Å" # lattice vectors are rows from cell matrix
    pbc = get(dict, "pbc", [true, true, true]) # default to periodic in all directions
    system_data[:boundary_conditions] =[p ? Periodic() : DirichletZero() for p in pbc]

    Atoms(NamedTuple(atom_data), NamedTuple(system_data))
end

read_dict(dict::Dict{String}{Any}) = Atoms(dict)

function write_dict(atoms::Atoms)
    system_data = Dict(pairs(atoms.system_data))
    atom_data = Dict(pairs(atoms.atom_data))
    if haskey(atom_data, :atomic_masses) 
        if atom_data[:atomic_masses] == getfield.(elements[atom_data[:atomic_numbers]], :atomic_mass)
            pop!(atom_data, :atomic_masses)
        elseif atom_data[:atomic_masses] == getfield.(elements[atom_data[:atomic_symbols]], :atomic_mass)
            pop!(atom_data, :atomic_masses)
        end
    end
    bcs = pop!(system_data, :boundary_conditions)
    box = pop!(system_data, :box)
    dict = Dict(
        "N_atoms" => length(atoms),
        "cell" => (_write_convert(box))',
        "pbc" => isequal.(bcs, [Periodic() for i=1:D]) |> Array,
        "info" => _dict_remap_rev(system_data),
        "arrays" => _dict_remap_rev(atom_data))
    return dict
end

# --------- AtomsBase interface

bounding_box(sys::Atoms)        = sys.system_data.box
boundary_conditions(sys::Atoms) = sys.system_data.boundary_conditions
Base.length(sys::Atoms)         = length(sys.atom_data.positions)
Base.size(sys::Atoms)           = size(sys.atom_data.positions)

function Base.isapprox(sys1::Atoms{NT1,NT2}, sys2::Atoms{NT1,NT2}) where {NT1, NT2}
    for (seq1, seq2) in [(sys1.system_data, sys2.system_data),
                         (sys1.atom_data, sys2.atom_data)]
        for (k1, v1) in pairs(seq1)
            if v1 isa Array{<:AbstractFloat} || v1 isa AbstractFloat
                v1 ≈ seq2[k1]  || (println("key $k1: $v1 !≈ $(seq2[k1])"); return false)
            else
                v1 == seq2[k1] || (println("key $k1: $v1 != $(seq2[k1])"); return false)
            end
        end
    end
    return true
end

# if types don't match then neither do the systems
Base.isapprox(sys1::Atoms{NT1,NT2}, sys2::Atoms{NT3,NT4}) where {NT1, NT2, NT3, NT4} = (println("type mismatch"); return false)

species_type(sys::FS) where {FS <: Atoms} = AtomView{FS}
Base.getindex(sys::Atoms, index::Int)     = AtomView(sys, index)

position(s::Atoms)       = s.atom_data.positions
atomic_symbol(s::Atoms)  = s.atom_data.atomic_symbols
atomic_number(s::Atoms)  = s.atom_data.atomic_numbers
atomic_mass(s::Atoms)    = s.atom_data.atomic_masses
velocity(s::Atoms)       = s.atom_data.velocities

position(s::Atoms, i)      = s.atom_data.positions[i]
atomic_symbol(s::Atoms, i) = s.atom_data.atomic_symbols[i]
atomic_number(s::Atoms, i) = s.atom_data.atomic_numbers[i]
atomic_mass(s::Atoms, i)   = s.atom_data.atomic_masses[i]
velocity(s::Atoms, i)      = s.atom_data.velocity[i]

# --------- FileIO compatible interface (hence not exported)

load(file::Union{String,IOStream}, frame; kwargs...) = Atoms(read_frame(file, frame; kwargs...))

function load(file::Union{String,IOStream}; kwargs...)
    seq = Atoms.(read_frames(file; kwargs...))
    if length(seq) == 1
        return seq[1]
    else
        return seq
    end
end

save(file::Union{String,IOStream}, system::Atoms; kwargs...) = write_frame(file, write_dict(system); kwargs...)

save(file::Union{String,IOStream}, systems::Vector{Atoms{NT1,NT2}}; kwargs...) where {NT1,NT2} = write_frames(file, write_dict.(systems); kwargs...)
