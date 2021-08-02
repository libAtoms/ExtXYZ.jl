using ExtXYZ
using Test

@testset "ExtXYZ.jl" begin
    @testset "read_extxyz" begin
        filename = "data/test.xyz"

        seq1 = ExtXYZ.read_frames(filename)
        seq2 = ExtXYZ.read_frames(filename, 4:10)
        frame = ExtXYZ.read_frames(filename, 4)
        @test all(seq1[4:10] .== seq2)
        
        f = open(filename, "r")
        seq3 = ExtXYZ.read_frames(f)
        close(f)
        
        @test all(seq1 .== seq3)    
    end
end
