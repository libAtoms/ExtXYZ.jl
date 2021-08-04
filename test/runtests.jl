using ExtXYZ
using Test

@testset "ExtXYZ.jl" begin
    @testset "read_extxyz" begin
        filename = "data/test.xyz"

        seq1 = read_frames(filename)
        seq2 = read_frames(filename, 4:10)
        @test all(seq1[4:10] .== seq2)

        frame4 = read_frame(filename, 4)
        @test frame4 == seq1[4]

        frame1 = read_frame(filename)
        @test frame1 == seq1[1]
        
        f = open(filename, "r")
        seq3 = read_frames(f)
        close(f)
        
        @test all(seq1 .== seq3)    
    end

    @testset "convert" begin
        filename = "data/test.xyz"
        nat, info, arrays = ExtXYZ.cfopen(filename, "r") do fp
            ExtXYZ.read_frame_dicts(fp)
        end
        # cinfo, values = convert(Ptr{ExtXYZ.DictEntry}, info)
    end
end
