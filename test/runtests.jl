using ExtXYZ
using JuLIP
using Test

@testset "ExtXYZ.jl" begin
    @testset "read_extxyz" begin
        filename = "data/test.xyz"

        seq1 = ExtXYZ.read_extxyz(filename)
        seq2 = ExtXYZ.read_extxyz(filename, 4:10)
        frame = ExtXYZ.read_extxyz(filename, 4)
        @test all(seq1[4:10] .== seq2)
        
        f = open(filename, "r")
        seq3 = ExtXYZ.read_extxyz(f)
        close(f)
        
        @test all(seq1 .== seq3)    
    end
end
