using ExtXYZ
using Test

# Memory-safety and performance sanity checks for the C interop layer
@testset "Stress" begin
    mktempdir() do dir
        path(name) = joinpath(dir, name)

        frame = Dict{String,Any}(
            "N_atoms" => 8,
            "cell" => [5.44 0.0 0.0; 0.0 5.44 0.0; 0.0 0.0 5.44],
            "pbc" => [true, true, true],
            "info" => Dict{String,Any}("energy" => -42.0, "config" => "bulk"),
            "arrays" => Dict{String,Any}(
                "species" => fill("Si", 8),
                "pos" => rand(3, 8)))

        @testset "repeated read/write with GC" begin
            # guards against double-free / use-after-free between Julia's
            # cfree_dict calls and the C library's internal error-path frees
            write_frame(path("stress.xyz"), frame)
            for i in 1:1000
                f = read_frame(path("stress.xyz"))
                @assert f["arrays"]["species"][1] == "Si"
                i % 100 == 0 && GC.gc()
            end
            # error path repeatedly (exercises free_partial_dicts in the C lib)
            write(path("bad.xyz"), "2\nLattice=\"oops Properties=species:S:1:pos:R:3\nSi 0 0 0\nGe 1 1 1\n")
            for i in 1:200
                @test_throws ErrorException read_frame(path("bad.xyz"))
                i % 50 == 0 && GC.gc()
            end
            @test true  # reached without crash
        end

        @testset "performance sanity" begin
            nframes = 1000
            frames = [frame for _ in 1:nframes]
            t_write = @elapsed write_frames(path("perf.xyz"), frames)
            t_read = @elapsed (seq = read_frames(path("perf.xyz")))
            t_read_tok = @elapsed read_frames(path("perf.xyz"); use_regex=false)
            @test length(seq) == nframes
            @info "performance sanity ($nframes frames, 8 atoms)" t_write t_read t_read_tok
            # generous bounds: just catch order-of-magnitude regressions
            @test t_write < 30
            @test t_read < 30
        end
    end
end
