using extxyz_jll

# On Windows, route fopen/fclose through libextxyz's wrappers so the FILE*
# is created and consumed by the same C runtime as extxyz_read_ll/extxyz_write_ll.
@static if Sys.iswindows()
    cfopen(filename::String, mode::String) = ccall((:extxyz_fopen, libextxyz),
                                                   Ptr{Cvoid},
                                                   (Cstring, Cstring),
                                                   filename, mode)

    function cfclose(fp::Ptr{Cvoid})
        (fp == C_NULL) && return
        ccall((:extxyz_fclose, libextxyz),
              Cint,
              (Ptr{Cvoid},),
              fp)
    end
else
    cfopen(filename::String, mode::String) = ccall(:fopen,
                                                   Ptr{Cvoid},
                                                   (Cstring, Cstring),
                                                   filename, mode)

    function cfclose(fp::Ptr{Cvoid})
        (fp == C_NULL) && return
        ccall(:fclose,
              Cint,
              (Ptr{Cvoid},),
              fp)
    end
end

function cfopen(f::Function, iostream::IOStream, mode::String="r")
    newfd = Libc.dup(RawFD(fd(iostream)))
    fp = ccall(@static(Sys.iswindows() ? :_fdopen : :fdopen),
               Ptr{Cvoid}, (Cint, Cstring), newfd, mode)
    try
        f(fp)
    finally
        cfclose(fp)
    end
end

function cfopen(f::Function, filename::String, mode::String="r")
    fp = cfopen(filename, mode)
    try
        f(fp)
    finally
        cfclose(fp)
    end
end

function cfopen(f::Function, iob::IOBuffer, mode::String="r")
    @static if Sys.iswindows()
        throw(ErrorException("Reading frames from IOBuffers is not supported on Windows."))
    end
    fp = ccall(:fmemopen, Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Cstring), pointer(iob.data, iob.ptr), iob.size, mode)
    try
        f(fp)
    finally
        cfclose(fp)
    end
end

struct DictEntry
    key::Ptr{Cchar}
    data::Ptr{Cvoid}
    data_t::Cint
    nrows::Cint
    ncols::Cint
    next::Ptr{DictEntry}
    first_data_ll::Ptr{Cvoid}
    last_data_ll::Ptr{Cvoid}
    n_in_row::Cint
end

const _kv_grammar = Ref{Ptr{Cvoid}}(0)

function __init__()
    _kv_grammar[] = ccall((:compile_extxyz_kv_grammar, libextxyz),
                           Ptr{Cvoid}, ())
    nothing
end

cfree_dict(dict::Ptr{DictEntry}) = ccall((:free_dict, libextxyz),
                                        Cvoid,
                                        (Ptr{DictEntry},),
                                        dict)

cprint_dict(dict::Ptr{DictEntry}) = ccall((:print_dict, libextxyz),
                                        Cvoid,
                                        (Ptr{DictEntry},),
                                        dict)


const DATA_I = 1
const DATA_F = 2
const DATA_B = 3
const DATA_S = 4

const TYPE_MAP = Dict(DATA_I => Cint,
                      DATA_F => Cdouble,
                      DATA_B => Cint,
                      DATA_S => Cstring)

import Base: convert

function convert(::Type{Dict{String,Any}}, c_dict::Ptr{DictEntry}; transpose_arrays=false)
    result = Dict{String,Any}()
    node_ptr = c_dict
    while node_ptr != C_NULL
        node = unsafe_load(node_ptr)
        data_ptr = reinterpret(Ptr{TYPE_MAP[node.data_t]}, node.data)

        if node.nrows == 0 && node.ncols == 0
            # scalar
            value = unsafe_load(data_ptr)
            # convert to primitive types
            if node.data_t == DATA_S
                value = unsafe_string(value)
            elseif node.data_t == DATA_B
                value = convert(Bool, value)
            end
        else
            # array, either 1D or 2D
            if node.nrows == 0
                # vector (1D array)
                dims = (node.ncols, )
            else
                # matrix (2D array)
                dims = (node.nrows, node.ncols)
            end

            if node.data_t == DATA_S && node.n_in_row < 0
                # libextxyz >= 0.4 returns per-atom string columns as a single
                # contiguous fixed-width NUL-padded buffer (cells in C row-major
                # order); n_in_row holds the negated cell width. A cell whose
                # string exactly fills the width has no NUL terminator, so the
                # search must be bounded by the cell - never unsafe_string().
                width = Int(-node.n_in_row)
                ncells = prod(dims)
                buf = unsafe_wrap(Array, Ptr{UInt8}(node.data), width * ncells)
                value = map(1:ncells) do i
                    cell = view(buf, (i-1)*width+1 : i*width)
                    nul = findfirst(iszero, cell)
                    String(cell[1:(nul === nothing ? width : nul-1)])
                end
                # same linear (C row-major) layout as the unsafe_wrap branch,
                # so the 2D reshape below applies unchanged
                node.nrows != 0 && (value = reshape(value, dims))
            else
                value = unsafe_wrap(Array,
                                    reinterpret(Ptr{TYPE_MAP[node.data_t]}, node.data),
                                    dims)

                if node.data_t == DATA_S
                    value = unsafe_string.(value)
                elseif node.data_t == DATA_B
                    value = convert(Array{Bool}, value)
                else
                    value = copy(value)
                end
            end

            if node.nrows != 0 && node.ncols != 0
                value = reshape(value, size(value, 2), size(value, 1))
                transpose_arrays && (value = permutedims(value, (2, 1)))
            end
        end

        key = unsafe_string(node.key)
        result[key] = value
        node_ptr = node.next
    end
    return result
end

# utility functions for converting from Julia types to corresponding C type and value

Ctype(::Type{S}) where S <: AbstractArray{T,N} where {T,N} = Ctype(T)
Ctype(::Type{<:Integer}) = (Cint, DATA_I)
Ctype(::Type{Bool}) = (Cint, DATA_B)
Ctype(::Type{<:Real}) = (Cdouble, DATA_F)
Ctype(::Type{String}) = (Cstring, DATA_S)

Cvalue(value::T; transpose_arrays=nothing) where {T<:Union{Integer,Bool,Real}} = convert(Ctype(T)[1], value)

function Cvalue(value::String; transpose_arrays=nothing) 
    ptr = pointer(Base.unsafe_convert(Cstring, Base.cconvert(Cstring, value)))
    new = Ptr{Cchar}(Libc.malloc(sizeof(value)+1))
    unsafe_copyto!(new, ptr, sizeof(value))
    unsafe_store!(new, C_NULL, sizeof(value)+1)
    return new
end

function Cvalue(value::Array{T,N}; transpose_arrays=false) where {T<:Union{Integer,Bool,Real},N}
    if ndims(value) == 2
        value = reshape(value, size(value, 2), size(value, 1))
    end
    value = convert(Array{Ctype(T)[1]}, value)
    if ndims(value) == 2 && transpose_arrays
        value = permutedims(value, (2, 1))
    end
    return value
end

function Cvalue(value::Array{String,1}; transpose_arrays=nothing)
    result = Array{Ptr{Cchar}}(undef, length(value))
    result .= Cvalue.(value)
    return result
end

# dims(value) -> (nrows, ncols) for scalar, vector and matrix value
dims(value) = (0, 0)
dims(value::AbstractVector) = (0, size(value, 1))
dims(value::AbstractMatrix) = (size(value, 1), size(value, 2))

function convert(::Type{Ptr{DictEntry}}, dict::Dict{String}{Any}; ordered_keys=nothing, transpose_arrays=false)
    c_dict_ptr = Ptr{DictEntry}(Libc.malloc(sizeof(DictEntry)))
    node_ptr = c_dict_ptr

    ordered_keys === nothing && (ordered_keys = keys(dict))
    for (idx, key) in enumerate(ordered_keys)
        value = dict[key]
        ckey = Cvalue(key)
        cvalue = Cvalue(value; transpose_arrays=transpose_arrays)
        nrow, ncol = dims(cvalue)
        type, data_t = Ctype(typeof(value))
        data = Ptr{type}(Libc.malloc(sizeof(cvalue)))
        if nrow == 0 && ncol == 0
            unsafe_store!(data, cvalue)
        else
            ptr = pointer(cvalue)
            data_t == DATA_S && (ptr = reinterpret(Ptr{type}, ptr)) # char** -> char*
            unsafe_copyto!(data, ptr, length(cvalue))
        end

        # allocate another DictEntry struct unless we're on the last one already
        next_ptr = C_NULL
        idx != length(dict) && (next_ptr = Ptr{DictEntry}(Libc.malloc(sizeof(DictEntry))))

        # n_in_row = 0 marks string data as the legacy char** layout (one
        # malloc per cell, as built by Cvalue above); the C writer and
        # free_dict only treat data as a contiguous buffer when n_in_row < 0
        node = DictEntry(ckey, data, data_t, nrow, ncol,
                         next_ptr, C_NULL, C_NULL, 0)
        unsafe_store!(node_ptr, node) # in place mutation of node_ptr
        node_ptr = next_ptr
    end
    return c_dict_ptr
end

# Error message emitted by libextxyz when the natoms header line cannot be
# parsed, quoting the offending line. Used to tell end-of-file apart from a
# genuine parse error: an empty buffer (true EOF) or a quoted line that is all
# whitespace (trailing blank lines) means EOF, anything else is a parse error.
# The text must match the C library verbatim (the Python bindings rely on the
# same message).
const _EOF_MESSAGE_RE = r"^Failed to parse int natoms from '(.*)'$"s

function _is_eof_message(msg::AbstractString)
    isempty(msg) && return true
    m = match(_EOF_MESSAGE_RE, msg)
    return m !== nothing && all(isspace, m[1])
end

function read_frame_dicts(fp::Ptr{Cvoid}; verbose=false, comment=nothing, use_regex=true)
    nat = Ref{Cint}(0)
    info = Ref{Ptr{DictEntry}}()
    arrays = Ref{Ptr{DictEntry}}()
    failed = false
    (comment === nothing) && (comment = C_NULL)
    # libextxyz writes diagnostics here with no NULL check - must be a real buffer
    error_message = zeros(UInt8, 1024)
    try
        res =  ccall((:extxyz_read_ll_opts, libextxyz),
                      Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ref{Cint}, Ptr{Ptr{DictEntry}}, Ptr{Ptr{DictEntry}}, Cstring, Ptr{UInt8}, Cint),
                      _kv_grammar[], fp, nat, info, arrays, comment, error_message, use_regex ? 0 : 1)
        if res != 1
            # on failure the C library frees the partial dicts itself, so the
            # finally block below must not free them again
            failed = true
            msg = GC.@preserve error_message unsafe_string(pointer(error_message))
            if _is_eof_message(msg)
                throw(EOFError())
            else
                error("extxyz parse error: $msg")
            end
        end
        if nat[] == 0
            failed = true
            error("ExtXYZ frame contains zero atoms. Behaviour is undefined!")
        end
        if verbose
            cprint_dict(info[])
            cprint_dict(arrays[])
        end

        pinfo = reinterpret(Ptr{DictEntry}, info[])
        parrays = reinterpret(Ptr{DictEntry}, arrays[])
        jinfo = convert(Dict{String,Any}, pinfo, transpose_arrays=true)
        jarrays = convert(Dict{String,Any}, parrays, transpose_arrays=false)
        return nat[], jinfo, jarrays

    finally
        if !failed
            cfree_dict(info[])
            cfree_dict(arrays[])
        end
    end
end

"""
extract "Lattice" entry and apply semantic conversions
"""
function extract_lattice!(result_dict)
    "Lattice" in keys(result_dict) || return nothing
    lattice = pop!(result_dict, "Lattice")
    if size(lattice) == (3, 3)
        lattice = convert(Array{Float64}, lattice)
    elseif size(lattice) == (3,)
        lattice = convert(Array{Float64}, diagm(lattice))
    elseif lattice.shape == (9,)
        lattice = convert(Array{Float64}, reshape(lattice, 3, 3))
    else
        error("Lattice has wrong shape!")
    end
    return lattice
end

"""
    read_frame(file; use_regex=true)

Read a single frame from the ExtXYZ file `file`, which can be a file pointer,
an open IO stream, a string filename or an IOBuffer.

Keyword arguments:
- `use_regex`: if `true` (default), per-atom lines are parsed with PCRE2 regular
  expressions. If `false`, a faster whitespace tokenizer is used instead; it
  still validates each field but is marginally more lenient on numeric formats.

Malformed input raises an `ErrorException` containing the parser's message.

Reading from IOBuffers is currently not supported on Windows.
"""
function read_frame(fp::Ptr{Cvoid}; verbose=false, kwargs...)
    nat, info, arrays = try
        read_frame_dicts(fp; verbose=verbose, kwargs...)
    catch err
        if isa(err, EOFError) 
            return nothing
        end
        rethrow()
    end

    dict = Dict{String, Any}()
    dict["N_atoms"] = nat # number of atoms

    # periodic boundary conditions
    if "pbc" in keys(info)
        dict["pbc"] = pop!(info, "pbc")
    end

    # cell is transpose of the stored lattice
    lattice = extract_lattice!(info)
    if (!isnothing(lattice)) dict["cell"] = permutedims(lattice, (2, 1)) end

    delete!(info, "Properties")

    # everything else stays in info and arrays
    dict["info"] = info
    dict["arrays"] = arrays

    return dict
end

read_frame(file::Union{String,IOStream,IOBuffer}, index; kwargs...) = only(read_frames(file, index; kwargs...))
read_frame(file::Union{String,IOStream,IOBuffer}; kwargs...) = read_frame(file, 1; kwargs...)

"""
    iread_frames(file[, range]; use_regex=true)

Return a Channel for reading from an ExtXYZ file. Frames are yielded one by one.
`range` can be a single integer, range object or integer array of frame indices.
See [`read_frame`](@ref) for the meaning of `use_regex`.

Example usage:

```julia
ch = iread_frames("file.xyz")
for frame in ch
    process(frame)
end
```
"""
function iread_frames(fp::Ptr{Cvoid}, range; close_fp=false, kwargs...)
    Channel() do channel
        for frame in 1:first(range)-1
            atoms = read_frame(fp; kwargs...)
            atoms === nothing && break
        end
        for frame in range
            atoms = read_frame(fp; kwargs...)
            atoms === nothing && break
            put!(channel, atoms)
        end
        if (close_fp)
            cfclose(fp)
        end
        channel
    end
end

function iread_frames(file::Union{String,IOStream,IOBuffer}, range; kwargs...)
    fp = cfopen(file, "r")
    fp == C_NULL && error("file $file cannot be opened for reading")
    iread_frames(fp, range; close_fp=true, kwargs...)
end

iread_frames(file::Union{String,IOStream,IOBuffer}, index::Int; kwargs...) = iread_frames(file, [index]; kwargs...)
iread_frames(file::Union{String,IOStream,IOBuffer}; kwargs...) = iread_frames(file, Iterators.countfrom(1); kwargs...)

"""
    read_frames(file[, range]; use_regex=true)

Read a sequence of frames from the ExtXYZ `file`, which can be specified by a file pointer, filename, IOStream or IOBuffer.

`range` can be a single integer, range object or integer array of frame indices.
See [`read_frame`](@ref) for the meaning of `use_regex`.

Reading from IOBuffers is currently not supported on Windows.
"""
read_frames(fp::Ptr{Cvoid}, range; kwargs...) = collect(iread_frames(fp, range; kwargs...))

function read_frames(file::Union{String,IOStream,IOBuffer}, range; kwargs...)
    cfopen(file) do fp
        fp == C_NULL && error("file $file cannot be opened for reading")
        read_frames(fp, range; kwargs...)
    end
end

read_frames(file::Union{String,IOStream,IOBuffer}, index::Int; kwargs...) = read_frames(file, [index]; kwargs...)
read_frames(file::Union{String,IOStream,IOBuffer}; kwargs...) = read_frames(file, Iterators.countfrom(1); kwargs...)

# C printf-style format string, or C_NULL to use the library default;
# strings are passed through so that ccall roots them for the call duration
_cfmt(fmt::AbstractString) = fmt
_cfmt(::Nothing) = Cstring(C_NULL)

function write_frame_dicts(fp::Ptr{Cvoid}, nat, info, arrays; verbose=false,
                           fmt_i=nothing, fmt_f=nothing, fmt_b=nothing, fmt_s=nothing)
    nat = Cint(nat)
    cinfo = convert(Ptr{DictEntry}, info; transpose_arrays=true)

    # ensure "species" (symbol!) goes in column 1 and "pos" goes in column 2
    ordered_keys = collect(keys(arrays))
    species_idx = findfirst(isequal("species"), ordered_keys)
    ordered_keys[1], ordered_keys[species_idx] = ordered_keys[species_idx], ordered_keys[1]
    pos_idx = findfirst(isequal("pos"), ordered_keys)
    ordered_keys[2], ordered_keys[pos_idx] = ordered_keys[pos_idx], ordered_keys[2]
    carrays = convert(Ptr{DictEntry}, arrays; ordered_keys=ordered_keys)

    if verbose
        cprint_dict(cinfo)
        cprint_dict(carrays)
    end
    try
        res = ccall((:extxyz_write_ll_fmt, libextxyz),
                     Cint, (Ptr{Cvoid}, Cint, Ptr{DictEntry}, Ptr{DictEntry}, Cstring, Cstring, Cstring, Cstring),
                     fp, nat, cinfo, carrays, _cfmt(fmt_i), _cfmt(fmt_f), _cfmt(fmt_b), _cfmt(fmt_s))
        res != 0 && error("error writing to file")
    finally
        cfree_dict(cinfo)
        cfree_dict(carrays)
    end
end

"""
    write_frame(file, dict; fmt_i=nothing, fmt_f=nothing, fmt_b=nothing, fmt_s=nothing)

Write a single atomic configuration represented by `dict` to `file`, which
can be a file pointer, open IO stream or string filename.

The `fmt_*` keyword arguments accept C printf-style format strings overriding
the output format for integer, float, bool and string values respectively
(e.g. `fmt_f="%21.16f"` for full double precision). When `nothing` (default),
the library defaults are used (`"%8d"`, `"%16.8f"`, `"%.1s"`, `"%s"`).
"""
function write_frame(fp::Ptr{Cvoid}, dict; kwargs...)
    nat = dict["N_atoms"]
    info = copy(dict["info"])
    if ("cell" in keys(dict)) info["Lattice"] = permutedims(dict["cell"], (2, 1)) end
    info["pbc"] = get(dict, "pbc", [true, true, true])

    write_frame_dicts(fp, nat, info, dict["arrays"]; kwargs...)
end

"""
    write_frames(file, dicts; append=false, fmt_i=nothing, fmt_f=nothing, fmt_b=nothing, fmt_s=nothing)

Write a sequence of atomic configurations to `file`. See [`write_frame`](@ref)
for the meaning of the `fmt_*` format-string keywords. Can also be used asynchronously
by passing a Channel in place of `dicts`, e.g.

```julia
Channel() do ch
    @async write_frames(outfile, ch)

    for frame in frames
        put!(ch, frame)
    end
end
```
"""
write_frames(fp::Ptr{Cvoid}, dicts; kwargs...) = write_frame.(dicts)

function write_frames(file::Union{String,IOStream}, dicts; append=false, kwargs...)
    mode = append ? "a" : "w"
    cfopen(file, mode) do fp
        fp == C_NULL && error("file $file cannot be opened for writing")
        for dict in dicts
            write_frame(fp, dict; kwargs...)
        end
    end
end

write_frame(file::Union{String,IOStream}, dict; kwargs...) = write_frames(file, [dict]; kwargs...)

export read_frame, read_frames, iread_frames, write_frame, write_frames
