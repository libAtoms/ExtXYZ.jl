using AtomsBase
using Unitful
using UnitfulAtomic
using StaticArrays

import AtomsBase: AbstractSystem
export Atoms

const D = 3 # TODO generalise to arbitrary spatial dimensions

# Types supported in the ExtXYZ C layer
const ExtxyzType = Union{Integer, AbstractFloat, AbstractString}

# ExtXYZ uses ASE units, see https://wiki.fysik.dtu.dk/ase/ase/units.html
# In particular note that uTime = u"Å" * sqrt(u"u" / u"eV") and thus
const uVelocity = sqrt(u"eV" / u"u")

"""
`struct Atoms` is the main type for flexible systems of atoms
"""
struct Atoms{P <: NamedTuple, Q <: NamedTuple} <: AbstractSystem{D}
    atom_data::P
    system_data::Q
end

function Atoms(system::AbstractSystem{D})
    n_atoms = length(system)
    s = species(system, :) 
    atomic_symbols = [Symbol(element(atomic_number(at)).symbol) for at in system]
    atomic_numbers = atomic_number.(s)
    if atomic_symbols != Symbol.(s)
        @warn("Mismatch between atomic numbers and atomic symbols, which is not supported " *
              "in ExtXYZ. Atomic numbers take preference.")
    end
    atom_data = Dict{Symbol,Any}(
        :atomic_symbol => atomic_symbols,
        :atomic_number => Int.(atomic_number(system, :)),  # gets messy if not Int
        :species => s, 
        :mass   => mass(system, :)
    )
    atom_data[:position] = map(1:n_atoms) do at
        pos = zeros(3)u"Å"
        pos[1:D] = position(system, at)
        SVector{D, eltype(pos)}(pos)  # AtomsBase 0.4 requires SVector
    end
    atom_data[:velocity] = map(1:n_atoms) do at
        vel = zeros(3) * uVelocity
        if !ismissing(velocity(system, :)) && !ismissing(velocity(system, at))
            vel[1:D] = velocity(system, at)
        end
        SVector{D, eltype(vel)}(vel)  # AtomsBase 0.4 requires SVector
    end

    for k in atomkeys(system)
        if k in (:species, :atomic_symbol, :atomic_number, :mass, :velocity, :position)
            continue  # Already done
        end
        # atomic_mass is deprecated for but is sometimes still used
        atoms_base_keys = (:charge, :atomic_mass, :covalent_radius, :vdw_radius,
                           :magnetic_moment, :pseudopotential)
        v = system[1, k]
        if k in atoms_base_keys || v isa ExtxyzType || v isa AbstractVector{<: ExtxyzType}
            # These are either Unitful quantities, which are uniformly supported
            # across all of AtomsBase or the value has a type that Extxyz can write
            # losslessly, so we can just write them no matter the value
            atom_data[k] = system[:, k]
        elseif v isa Quantity || (v isa AbstractArray && eltype(v) <: Quantity)
            @warn "Unitful quantity $k is not yet supported in ExtXYZ."
        else
            @warn "Writing quantities of type $(typeof(v)) is not supported in ExtXYZ."
        end
    end

    system_data = Dict{Symbol,Any}(
        :bounding_box => bounding_box(system),
        :periodicity => periodicity(system) 
    )

    # Extract extra system properties
    for (k, v) in pairs(system)
        atoms_base_keys = (:charge, :multiplicity, :periodicity, :bounding_box)
        if k in atoms_base_keys || v isa ExtxyzType || v isa AbstractArray{<: ExtxyzType}
            # These are either Unitful quantities, which are uniformly supported
            # across all of AtomsBase or the value has a type that Extxyz can write
            # losslessly, so we can just write them no matter the value
            system_data[k] = v
        elseif v isa Quantity || (v isa AbstractArray && eltype(v) <: Quantity)
            @warn "Unitful quantity $k is not yet supported in ExtXYZ."
        else
            @warn "Writing quantities of type $(typeof(v)) is not supported in ExtXYZ."
        end
    end

    ExtXYZ.Atoms(NamedTuple(atom_data), NamedTuple(system_data))
end


function Atoms(dict::Dict{String, Any})
    arrays = dict["arrays"]
    info   = dict["info"]

    if haskey(arrays, "Z")
        Z = Int.(arrays["Z"])
    elseif haskey(arrays, "species")
        Z = [element(Symbol(spec)).number for spec in arrays["species"]]
    else
        error("Cannot determine atomic numbers. Either 'Z' or 'species' must " *
              "be present in arrays")
    end
    @assert length(Z) == dict["N_atoms"]

    atom_data = Dict{Symbol, Any}(
        :position      => collect(eachcol(arrays["pos"]))u"Å",
        :atomic_number => Z,
    )
    if haskey(arrays, "species")
        atom_data[:atomic_symbol] = Symbol.(arrays["species"])
    else
        atom_data[:atomic_symbol] = [Symbol(element(num).symbol) for num in Z]
    end
    if haskey(arrays, "mass")
        atom_data[:mass] = arrays["mass"]u"u"
    else
        atom_data[:mass] = [element(num).atomic_mass for num in Z]
    end
    if haskey(arrays, "velocities")
        atom_data[:velocity] = collect(eachcol(arrays["velocities"])) * uVelocity
    else
        atom_data[:velocity] = [zeros(3) * uVelocity for _ in Z]
    end

    for key in keys(arrays)
        key in ("mass", "species", "Z", "pos", "velocities") && continue  # Already done
        if key in ("vdw_radius", "covalent_radius")  # Add length unit
            atom_data[Symbol(key)] = arrays[key] * u"Å"
        elseif key in ("charge", )  # Add charge unit
            atom_data[Symbol(key)] = arrays[key] * u"e_au"
        elseif typeof(arrays[key]) <: AbstractMatrix
            atom_data[Symbol(key)] = [ collect(col) for col in eachcol(arrays[key]) ]
        else
            atom_data[Symbol(key)] = arrays[key]
        end
    end

    system_data = Dict{Symbol, Any}()
    if haskey(dict, "cell")
        system_data[:bounding_box] = collect(eachrow(dict["cell"]))u"Å"
        if haskey(dict, "pbc")
            system_data[:periodicity] = tuple(dict["pbc"]...)
        else
            @warn "'pbc' not contained in dict. Defaulting to all-periodic boundary. "
            system_data[:periodicity] = (true, true, true)
        end
    else  # Infinite system
        haskey(dict, "pbc") && @warn "'pbc' ignored since no 'cell' entry found in dict."
        system_data[:periodicity] = (false, false, false) 
        system_data[:bounding_box] = ( SVector(Inf, 0.0, 0.0), 
                                       SVector(0.0, Inf, 0.0), 
                                       SVector(0.0, 0.0, Inf) )
    end

    for key in keys(info)
        if key in ("charge", )
            system_data[Symbol(key)] = info[key] * u"e_au"
        else
            system_data[Symbol(key)] = info[key]
        end
    end

    Atoms(NamedTuple(atom_data), NamedTuple(system_data))
end

read_dict(dict::Dict{String,Any}) = Atoms(dict)

function write_dict(atoms::Atoms)
    arrays = Dict{String,Any}()

    arrays["Z"] = atoms.atom_data.atomic_number
    arrays["atomic_symbol"] = [element(Z).symbol for Z in arrays["Z"]]
    if atoms.atom_data.atomic_symbol != [Symbol(element(Z).symbol) for Z in arrays["Z"]]
        @warn("Mismatch between atomic numbers and atomic symbols, which is not supported " *
              "in ExtXYZ. Atomic numbers take preference.")
    end
    if atoms.atom_data.mass != [element(Z).atomic_mass for Z in arrays["Z"]]
        arrays["mass"] = ustrip.(u"u", atoms.atom_data.mass)
    end

    arrays["velocities"] = zeros(D, length(atoms))
    for (i, velocity) in enumerate(atoms.atom_data.velocity)
        arrays["velocities"][:, i] = ustrip.(uVelocity, velocity)
    end
    arrays["pos"] = zeros(D, length(atoms))
    for (i, position) in enumerate(atoms.atom_data.position)
        arrays["pos"][:, i] = ustrip.(u"Å", position)
    end

    for (k, v) in pairs(atoms.atom_data)
        k in (:mass, :atomic_mass, :atomic_symbol, :atomic_number, :position, 
              :velocity, :species) && continue
        if k in (:vdw_radius, :covalent_radius)  # Remove length unit
            arrays[string(k)] = ustrip.(u"Å", v)
        elseif k in (:charge, )
            arrays[string(k)] = ustrip.(u"e_au", v)
        elseif v isa AbstractVector{<:ExtxyzType}
            arrays[string(k)] = v  # These can be written losslessly
        elseif v isa AbstractArray && eltype(v) <: AbstractVector{<:ExtxyzType}
            arrays[string(k)] = reduce(hcat, v)
        else
            @warn "Writing quantities of type $(typeof(v)) is not supported in write_dict."
        end
    end

    pbc = atoms.system_data.periodicity
    cell = zeros(D, D)
    for (i, bvector) in enumerate(atoms.system_data.bounding_box)
        cell[i, :] = ustrip.(u"Å", bvector)
    end

    # Deal with other system keys
    info = Dict{String,Any}()
    for (k, v) in pairs(atoms.system_data)
        k in (:periodicity, :bounding_box) && continue # Already dealt with
        if k in (:charge, )
            info[string(k)] = ustrip(u"e_au", atoms.system_data[k])
        elseif v isa ExtxyzType
            info[string(k)] = v
        elseif v isa AbstractArray{<:ExtxyzType}
            if size(v, 1) != 3
                @warn("Writing Array data with a first leading dimension " *
                      "different from 3 not supported.")  # Else leads to segfaults
            else
                info[string(k)] = convert(Array, v)
            end
        else
            @warn "Writing quantities of type $(typeof(v)) is not supported in write_dict."
        end
    end
    dict = Dict("N_atoms" => length(atoms),
                "pbc"     => [pbc...],
                "info"    => info,
                "arrays"  => arrays)
    all(cell .!= Inf) && (dict["cell"] = cell) # only write cell if its finite
    return dict
end
write_dict(system::AbstractSystem{D}) = write_dict(Atoms(system))

# --------- AtomsBase interface

Base.length(sys::Atoms) = length(sys.atom_data.position)
Base.size(sys::Atoms)   = (length(sys), )
AtomsBase.bounding_box(sys::Atoms) = sys.system_data.bounding_box
AtomsBase.periodicity(sys::Atoms) = sys.system_data.periodicity

# AtomsBase now requires a cell object to be returned instead of bounding_box 
# and boundary conditions. But this can just be constructed on the fly. 
AtomsBase.cell(sys::Atoms) = AtomsBase.PeriodicCell(; 
                cell_vectors = sys.system_data.bounding_box, 
                 periodicity = sys.system_data.periodicity )

Base.getindex(sys::Atoms, x::Symbol) = getindex(sys.system_data, x)
Base.haskey(sys::Atoms, x::Symbol)   = haskey(sys.system_data, x)
Base.keys(sys::Atoms) = keys(sys.system_data)

Base.getindex(sys::Atoms, i::Integer) = AtomView(sys, i)
Base.getindex(sys::Atoms, i::Integer, x::Symbol) = getindex(sys.atom_data, x)[i]
Base.getindex(sys::Atoms, i::AbstractVector{<: Integer}, x::Symbol) = getindex(sys.atom_data, x)[i]
Base.getindex(sys::Atoms, ::Colon,    x::Symbol) = getindex(sys.atom_data, x)

AtomsBase.atomkeys(sys::Atoms) = keys(sys.atom_data)
AtomsBase.hasatomkey(sys::Atoms, x::Symbol) = haskey(sys.atom_data, x)

const _IDX = Union{Colon, Integer, AbstractArray{<: Integer}}
AtomsBase.position(s::Atoms, i::_IDX)      = getindex(s, i, :position)
AtomsBase.velocity(s::Atoms, i::_IDX)      = getindex(s, i, :velocity)
AtomsBase.mass(s::Atoms, i::_IDX)          = getindex(s, i, :mass)
AtomsBase.atomic_symbol(s::Atoms, i::_IDX) = getindex(s, i, :atomic_symbol)
AtomsBase.atomic_number(s::Atoms, i::_IDX) = getindex(s, i, :atomic_number)

# AtomsBase now requires the `species` function to be implemented. Since 
# ExtXYZ requires that atoms are uniquely identified by their atomic number, we 
# will use the atomic number as the species identifier.
AtomsBase.species(s::Atoms, i::_IDX) = 
            AtomsBase.ChemicalSpecies.(AtomsBase.atomic_symbol(s, i))

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

save(file::Union{String,IOStream}, system::AbstractSystem; kwargs...) = write_frame(file, write_dict(system); kwargs...)
save(file::Union{String,IOStream}, systems::AbstractVector{<: AbstractSystem}; kwargs...) = write_frames(file, write_dict.(systems); kwargs...)
