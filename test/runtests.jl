using ExtXYZ
using Test
using PyCall
using Unitful
using AtomsBase

"""
Check two dictionaries `seq1` and `seq2` are approximately equal

Returns true if all keys match, all non-float values match, and all pairs 
of float values `(v1, v2)` satisfy `isapprox(v1, v2)`
"""
function Base.isapprox(seq1::AbstractDict, seq2::AbstractDict)
    for (k1,v1) in seq1
        k1 ∈ keys(seq2) || (println("key $k1 missing from seq2"); return false)
        if v1 isa AbstractDict
            v1 ≈ seq2[k1] || (println("key $k1: $v1 !≈ $(seq2[k1])"); return false)
        elseif v1 isa Array{<:AbstractFloat} || v1 isa AbstractFloat
            v1 ≈ seq2[k1]  || (println("key $k1: $v1 !≈ $(seq2[k1])"); return false)
        else
            v1 == seq2[k1] || (println("key $k1: $v1 != $(seq2[k1])"); return false)
        end
    end
    return true
end

function Test_Units(seq::ExtXYZ.Atoms)
    all([all(unit.(position(seq)[i]) .== u"Å") for i=1:length(seq)]) || return false
    all([all(unit.(bounding_box(seq)[i]) .== u"Å") for i=1:n_dimensions(seq)]) || return false
    all([all(unit.(atomic_mass(seq)[i]) .== u"u") for i=1:length(seq)]) || return false
    return true
end

@testset "ExtXYZ.jl" begin
    test_xyz(frame) = """2
cutoff=5.50000000 pbc=[T, T, T] nneightol=1.20000000 config_type=isolated_atom matrix=[[$frame.00000000, 1.00000000, 2.00000000], [3.00000000, 4.00000000, 5.00000000], [6.00000000, 7.00000000, 8.00000000]] dft_energy=-158.$(frame)4496821 Lattice="1.00000000 2.00000000 3.00000000 4.00000000 5.00000000 6.00000000 7.00000000 $(frame).00000000 9.00000000" gap_energy=-157.72725320  Properties=species:S:1:pos:R:3:n_neighb:I:1:dft_force:R:3:gap_force:R:3:map_shift:I:3
Si        10.00000000      11.00000000      $frame.00000000          0         0.10000000       0.20000000       0.30000000         0.1$(frame)000000       0.22000000       0.330000000         -1       -2       0
Si        13.00000000      14.00000000      $(frame+1).00000000          0         0.10000000       0.20000000       0.30000000         0.2$(frame)000000       0.22000000       0.330000000         -1       -2       0
"""

    infile = "test.xyz"
    open(infile, "w") do io
        for frame=1:10
            print(io, test_xyz(frame))
        end
    end
    outfile = "dump.xyz"

    try
        @testset "read" begin
            seq1 = read_frames(infile)

            @test all([s["cell"][3,2] for s in seq1] .== 1.0:10.0)
            @test all([s["arrays"]["pos"][3,1] for s in seq1] .== 1.0:10.0)
            @test all([s["arrays"]["pos"][3,2] for s in seq1] .== 2.0:11.0)
            @test all([s["info"]["matrix"][1,1] for s in seq1] .== 1.0:10.0)
            @test all([s["info"]["matrix"][2,1] for s in seq1] .== 3.0)

            @test all([s["info"]["dft_energy"] for s in seq1] .== [parse(Float64, "-158.$(frame)4496821") for frame=1:10])
            @test all([s["arrays"]["gap_force"][1,1] for s in seq1] .== [parse(Float64, "0.1$(frame)") for frame=1:10])
            @test all([s["arrays"]["gap_force"][1,2] for s in seq1] .== [parse(Float64, "0.2$(frame)") for frame=1:10])

            seq2 = read_frames(infile, 4:10)
            @test all(seq1[4:10] .== seq2)

            frame4 = read_frame(infile, 4)
            @test frame4 == seq1[4]

            frame1 = read_frame(infile)
            @test frame1 == seq1[1]
            
            f = open(infile, "r")
            seq3 = read_frames(f)
            close(f)
            
            @test all(seq1 .== seq3)
        end

        @testset "convert" begin
            nat, info, arrays = ExtXYZ.cfopen(infile, "r") do fp
                ExtXYZ.read_frame_dicts(fp)
            end
            cinfo = convert(Ptr{ExtXYZ.DictEntry}, info)
            carrays = convert(Ptr{ExtXYZ.DictEntry}, arrays)

            ninfo = convert(Dict{String,Any}, cinfo)
            narrays = convert(Dict{String,Any}, carrays)

            @test ninfo == info
            @test narrays == arrays
        end

        @testset "write" begin
            seq1 = read_frames(infile)
            write_frames(outfile, seq1)
            seq2 = read_frames(outfile)
            @test all(seq1 .≈ seq2) # use custom isapprox(), since expect some loss of precision on round-trip
        end

        @testset "iread" begin
            seq1 = read_frames(infile)
            seq2 = collect(iread_frames(infile))
            @test all(seq1 .≈ seq2)
        end

        @testset "iwrite" begin
            seq1 = read_frames(infile)
            ch = Channel()
            job = @async write_frames(outfile, ch)
            for frame in seq1
                put!(ch, frame)
            end
            close(ch)
            wait(job)
            seq2 = read_frames(outfile)
            @test all(seq1 .≈ seq2)
        end

        @testset "AtomsBase" begin
            seq1 = ExtXYZ.load(infile)
            ExtXYZ.save(outfile, seq1)
            seq2 = ExtXYZ.load(outfile)
            @test all(seq1 .≈ seq2)
            @test all(repr(seq1) == repr(seq2))

            @testset "AtomsBase unit_tests" begin
                at = ExtXYZ.load(infile, 1)
                @test at ≈ seq1[1]
                @test all(Test_Units.(vcat(seq1, seq2)))
                atom_data = Dict(:positions => [[1,2,3],[2,3,4]].*10^(-12)*u"m", :velocities => [[3,4,5],[4,5,6]]*u"Å/(Å*sqrt(u/eV))", :atomic_numbers => [25, 25], :atomic_symbols => [:Mn, :Mn])
                system_data = Dict(:box => [[3,3,3],[4,4,4],[5,5,5]]*u"Å", :boundary_conditions => [Periodic(), Periodic(), Periodic()], 
                :dft_energy => 10*u"eV", :Momentum => 11*u"u*Å/(Å*sqrt(u/eV))", :Force => 12*u"eV/Å", :Stress => [13,14]*u"eV/Å^3", :conv_vec => [10,11,12].*10^(-12)*u"m", :conv_val => 13*10^(-12)*u"m")
                sys1 = Atoms(NamedTuple(atom_data), NamedTuple(system_data))
                ExtXYZ.save(outfile, sys1)
                sys2 = ExtXYZ.load(outfile)
                @test(position(sys2) == [[0.010,0.020,0.030],[0.020,0.030,0.040]]*u"Å")
                @test(velocity(sys2) == [[3,4,5],[4,5,6]]*u"Å/(Å*sqrt(u/eV))")
                @test(sys2.system_data[:dft_energy] == 10)
                @test(sys2.system_data[:Momentum] == 11)
                @test(sys2.system_data[:Force] == 12)
                @test(sys2.system_data[:Stress] == [13,14])
                @test(sys2.system_data[:conv_vec] == [0.10, 0.11, 0.12])
                @test(sys2.system_data[:conv_val] == 0.13)
            end
        end

        try
            ase_io = pyimport("ase.io")

            @testset "ASE" begin
                seq = read_frames(infile)
                ase_seq = ase_io.read(infile * "@:")
    
                for (frame, ase_atoms) in zip(seq, ase_seq)
                    @test frame["N_atoms"] ≈ length(ase_atoms)
                    @test frame["arrays"]["pos"] ≈ ase_atoms.positions'
                    @test frame["arrays"]["gap_force"] ≈ ase_atoms.arrays["gap_force"]'
                    @test frame["arrays"]["n_neighb"] ≈ ase_atoms.arrays["n_neighb"]
                    @test frame["cell"] ≈ ase_atoms.cell.array
                end
            end            
        catch e
            if e isa PyError
                println("ASE not installed, skipping ASE comparison tests")
            else
                throw(e)
            end
        end

    finally
        rm(infile, force=true)
        rm(outfile, force=true)
    end
end
