using AtomsBase
using Unitful
using UnitfulAtomic

export Atoms

const D = 3 # TODO generalise to arbitrary spatial dimensions

# Types supported in the ExtXYZ C layer
const ExtxyzType = Union{Integer, AbstractFloat, AbstractString}

"""
`struct Atoms` is the main type for flexible systems of atoms
"""
struct Atoms{P <: NamedTuple, Q <: NamedTuple} <: AbstractSystem{D}
    atom_data::P
    system_data::Q
end

function Atoms(system::AbstractSystem{D})
    n_atoms = length(system)
    atomic_symbols = [Symbol(element(atomic_number(at)).symbol) for at in system]
    if atomic_symbols != atomic_symbol(system)
        @warn("Mismatch between atomic numbers and atomic symbols, which is not supported " *
              "in ExtXYZ. Atomic numbers take preference.")
    end
    atom_data = Dict{Symbol,Any}(
        :atomic_symbol => atomic_symbols,
        :atomic_number => atomic_number(system),
        :atomic_mass   => atomic_mass(system)
    )
    atom_data[:position] = map(1:n_atoms) do at
        pos = zeros(3)u"Å"
        pos[1:D] = position(system, at)
        pos
    end
    atom_data[:velocity] = map(1:n_atoms) do at
        vel = zeros(3)u"Å/s"
        if !ismissing(velocity(system)) && !ismissing(velocity(system, at))
            vel[1:D] = velocity(system, at)
        end
        vel
    end

    for k in atomkeys(system)
        if k in (:atomic_symbol, :atomic_number, :atomic_mass, :velocity, :position)
            continue  # Already done
        end
        atoms_base_keys = (:charge, :covalent_radius, :vdw_radius,
                           :magnetic_moment, :pseudopotential)
        v = system[1, k]
        if k in atoms_base_keys || v isa ExtxyzType
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

    box = map(1:3) do i
        v = zeros(3)u"Å"
        i ≤ D && (v[1:D] = bounding_box(system)[i])
        v
    end
    system_data = Dict{Symbol,Any}(
        :bounding_box => box,
        :boundary_conditions => boundary_conditions(system)
    )

    # Extract extra system properties
    system_data = Dict{Symbol,Any}()
    for (k, v) in pairs(system)
        atoms_base_keys = (:charge, :multiplicity, :boundary_conditions, :bounding_box)
        if k in atoms_base_keys || v isa ExtxyzType
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
        atom_data[:atomic_mass] = arrays["mass"]u"u"
    else
        atom_data[:atomic_mass] = [element(num).atomic_mass for num in Z]
    end
    if haskey(arrays, "velocities")
        atom_data[:velocity] = collect(eachcol(arrays["velocities"]))u"Å/s"
    else
        atom_data[:velocity] = [zeros(3)u"Å/s" for _ in 1:Z]
    end

    for key in keys(arrays)
        key in ("mass", "species", "Z", "pos", "velocities") && continue  # Already done
        if key in ("vdw_radius", "covalent_radius")  # Add length unit
            atom_data[Symbol(key)] = arrays[key] * u"Å"
        elseif key in ("charge", )  # Add charge unit
            atom_data[Symbol(key)] = arrays[key] * u"e_au"
        else
            atom_data[Symbol(key)] = arrays[key]
        end
    end

    system_data = Dict{Symbol, Any}(:bounding_box => collect(eachrow(dict["cell"]))u"Å", )
    if haskey(dict, "pbc")
        system_data[:boundary_conditions] = [p ? Periodic() : DirichletZero()
                                             for p in dict["pbc"]]
    else
        @warn "'pbc' not contained in 'info' dict. Defaulting to all-periodic boundary. "
        system_data[:boundary_conditions] = fill(Periodic(), 3)
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
    arrays["species"] = [element(Z).symbol for Z in arrays["Z"]]
    if atoms.atom_data.atomic_symbol != [Symbol(element(Z).symbol) for Z in arrays["Z"]]
        @warn("Mismatch between atomic numbers and atomic symbols, which is not supported " *
              "in ExtXYZ. Atomic numbers take preference.")
    end
    if atoms.atom_data.atomic_mass != [element(Z).atomic_mass for Z in arrays["Z"]]
        arrays["mass"] = ustrip.(u"u", atoms.atom_data.atomic_mass)
    end

    arrays["velocities"] = zeros(D, length(atoms))
    for (i, velocity) in enumerate(atoms.atom_data.velocity)
        arrays["velocities"][:, i] = ustrip.(u"Å/s", velocity)
    end
    arrays["pos"] = zeros(D, length(atoms))
    for (i, position) in enumerate(atoms.atom_data.position)
        arrays["pos"][:, i] = ustrip.(u"Å", position)
    end

    for (k, v) in pairs(atoms.atom_data)
        k in (:atomic_mass, :atomic_symbol, :atomic_number, :position, :velocity) && continue
        if k in (:vdw_radius, :covalent_radius)  # Remove length unit
            arrays[string(k)] = ustrip.(u"Å", v)
        elseif k in (:charge, )
            arrays[string(k)] = ustrip.(u"e_au", v)
        elseif v isa AbstractVector{<:ExtxyzType}
            arrays[string(k)] = v  # These can be written losslessly
        else
            @warn "Writing quantities of type $(typeof(v)) is not supported in write_dict."
        end
    end

    pbc = zeros(Bool, D)
    for (i, bc) in enumerate(atoms.system_data.boundary_conditions)
        pbc[i] = bc isa Periodic
    end
    cell = zeros(D, D)
    for (i, bvector) in enumerate(atoms.system_data.bounding_box)
        cell[i, :] = ustrip.(u"Å", bvector)
    end

    # Deal with other system keys
    info = Dict{String,Any}()
    for (k, v) in pairs(atoms.system_data)
        k in (:boundary_conditions, :bounding_box) && continue # Already dealt with
        if k in (:charge, )
            info[string(k)] = ustrip(u"e_au", atoms.system_data[k])
        elseif v isa ExtxyzType
            info[string(k)] = v
        elseif v isa AbstractArray{<:ExtxyzType}
            info[string(k)] = convert(Array, v)
        else
            @warn "Writing quantities of type $(typeof(v)) is not supported in write_dict."
        end
    end
    Dict("N_atoms" => length(atoms),
         "cell"    => cell,
         "pbc"     => pbc,
         "info"    => info,
         "arrays"  => arrays)
end
write_dict(system::AbstractSystem{D}) = write_dict(Atoms(system))

# --------- AtomsBase interface

Base.length(sys::Atoms) = length(sys.atom_data.position)
Base.size(sys::Atoms)   = size(sys.atom_data.position)
AtomsBase.bounding_box(sys::Atoms) = sys.system_data.bounding_box
AtomsBase.boundary_conditions(sys::Atoms) = sys.system_data.boundary_conditions

AtomsBase.species_type(::FS) where {FS <: Atoms} = AtomView{FS}
Base.getindex(sys::Atoms, x::Symbol) = getindex(sys.system_data, x)
Base.haskey(sys::Atoms, x::Symbol)   = haskey(sys.system_data, x)
Base.keys(sys::Atoms) = keys(sys.system_data)

Base.getindex(sys::Atoms, i::Integer) = AtomView(sys, i)
Base.getindex(sys::Atoms, i::Integer, x::Symbol) = getindex(sys.atom_data, x)[i]
Base.getindex(sys::Atoms, ::Colon,    x::Symbol) = getindex(sys.atom_data, x)

AtomsBase.atomkeys(sys::Atoms) = keys(sys.atom_data)
AtomsBase.hasatomkey(sys::Atoms, x::Symbol) = haskey(sys.atom_data, x)

AtomsBase.position(s::Atoms)                  = Base.getindex(s, :, :position)
AtomsBase.position(s::Atoms, i::Integer)      = Base.getindex(s, i, :position)
AtomsBase.velocity(s::Atoms)                  = Base.getindex(s, :, :velocity)
AtomsBase.velocity(s::Atoms, i::Integer)      = Base.getindex(s, i, :velocity)
AtomsBase.atomic_mass(s::Atoms)               = Base.getindex(s, :, :atomic_mass)
AtomsBase.atomic_mass(s::Atoms, i::Integer)   = Base.getindex(s, i, :atomic_mass)
AtomsBase.atomic_symbol(s::Atoms)             = Base.getindex(s, :, :atomic_symbol)
AtomsBase.atomic_symbol(s::Atoms, i::Integer) = Base.getindex(s, i, :atomic_symbol)
AtomsBase.atomic_number(s::Atoms)             = Base.getindex(s, :, :atomic_number)
AtomsBase.atomic_number(s::Atoms, i::Integer) = Base.getindex(s, i, :atomic_number)

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
