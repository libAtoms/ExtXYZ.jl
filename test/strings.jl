using ExtXYZ
using Test

# Per-atom string columns: libextxyz >= 0.4 returns these as a contiguous
# fixed-width buffer (initial cell width 8, grown on demand), unlike the
# legacy char** layout still used for info-dict strings and Julia-built dicts.
@testset "String columns" begin
    mktempdir() do dir
        path(name) = joinpath(dir, name)

        @testset "width growth mid-column" begin
            # the long name appears LATE in the frame, forcing the C parser to
            # reallocate the column buffer after earlier cells are filled
            species = ["H", "C", "Si", "Cl", "VeryLongSpeciesName", "O"]
            lines = join(["$s 0.0 0.0 $(i-1).0" for (i, s) in enumerate(species)], "\n")
            write(path("long.xyz"), """$(length(species))
Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44" Properties=species:S:1:pos:R:3
$lines
""")
            for use_regex in (true, false)
                frame = read_frame(path("long.xyz"); use_regex=use_regex)
                @test frame["arrays"]["species"] == species
                @test frame["arrays"]["pos"][3, :] == 0.0:(length(species)-1)
            end
        end

        @testset "cell boundary lengths" begin
            # lengths around the initial cell width of 8: 7 fills a cell up to
            # its final NUL, 8 forces growth to the next multiple of 8
            species = ["Abcdefg", "Abcdefgh", "Abcdefghijklmno", "H"]
            lines = join(["$s 0.0 0.0 0.0" for s in species], "\n")
            write(path("widths.xyz"), """$(length(species))
Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44" Properties=species:S:1:pos:R:3
$lines
""")
            frame = read_frame(path("widths.xyz"))
            @test frame["arrays"]["species"] == species
        end

        @testset "multi-column string property" begin
            # tags:S:2 -> per-atom string matrix; orientation must match the
            # numeric convention (N_component x N_atoms, like pos being 3 x nat)
            write(path("tags.xyz"), """3
Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44" Properties=species:S:1:pos:R:3:tags:S:2
Si 0.0 0.0 0.0 a1 b1
Ge 1.0 1.0 1.0 a2 b2
C  2.0 2.0 2.0 a3 b3
""")
            frame = read_frame(path("tags.xyz"))
            @test frame["arrays"]["tags"] == ["a1" "a2" "a3";
                                              "b1" "b2" "b3"]
        end

        @testset "info-dict string array (legacy layout)" begin
            write(path("infostr.xyz"), """1
Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44" names=[alpha, beta, gamma] Properties=species:S:1:pos:R:3
Si 0.0 0.0 0.0
""")
            frame = read_frame(path("infostr.xyz"))
            @test frame["info"]["names"] == ["alpha", "beta", "gamma"]
        end

        @testset "round trips" begin
            # C-read -> Julia-write -> C-read idempotence with long species
            # (the writer adds a default "pbc" key, so compare frame's keys
            # against the round-tripped result, not the other way around)
            frame = read_frame(path("long.xyz"))
            write_frame(path("rt.xyz"), frame)
            @test frame ≈ read_frame(path("rt.xyz"))

            # Julia-built frame (legacy char** write path) -> C-read
            nat = 3
            built = Dict{String,Any}(
                "N_atoms" => nat,
                "cell" => [5.44 0.0 0.0; 0.0 5.44 0.0; 0.0 0.0 5.44],
                "pbc" => [true, true, true],
                "info" => Dict{String,Any}("comment" => "test"),
                "arrays" => Dict{String,Any}(
                    "species" => ["H", "LongSpeciesName", "O"],
                    "pos" => [0.0 1.0 2.0; 0.0 1.0 2.0; 0.0 1.0 2.0]))
            write_frame(path("built.xyz"), built)
            frame2 = read_frame(path("built.xyz"))
            @test frame2["arrays"]["species"] == built["arrays"]["species"]
            @test frame2["arrays"]["pos"] ≈ built["arrays"]["pos"]
        end

        @testset "tokenizer equivalence" begin
            frames = [Dict{String,Any}(
                          "N_atoms" => 2,
                          "cell" => [5.44 0.0 0.0; 0.0 5.44 0.0; 0.0 0.0 5.44],
                          "pbc" => [true, true, true],
                          "info" => Dict{String,Any}("step" => i, "energy" => -1.5i),
                          "arrays" => Dict{String,Any}(
                              "species" => ["Si", "Ge"],
                              "pos" => [0.0 1.0; 0.0 1.0; 0.0 Float64(i)]))
                      for i in 1:5]
            write_frames(path("traj.xyz"), frames)
            seq_regex = read_frames(path("traj.xyz"); use_regex=true)
            seq_token = read_frames(path("traj.xyz"); use_regex=false)
            @test seq_regex == seq_token
        end

        @testset "comment-line parser equivalence" begin
            # the first-char-dispatch comment parser (use_cleri=false) must be
            # bit-identical to the libcleri grammar over assorted info types
            frames = [Dict{String,Any}(
                          "N_atoms" => 2,
                          "cell" => [5.44 0.0 0.0; 0.0 5.44 0.0; 0.0 0.0 5.44],
                          "pbc" => [true, true, false],
                          "info" => Dict{String,Any}("step" => i, "energy" => -1.5i,
                                                     "label" => "frame $i", "ok" => isodd(i),
                                                     "vec" => Float64[i, 2i, 3i]),
                          "arrays" => Dict{String,Any}(
                              "species" => ["Si", "O"],
                              "pos" => [0.0 1.0; 0.0 1.0; 0.0 Float64(i)]))
                      for i in 1:5]
            write_frames(path("ct.xyz"), frames)
            # all four (use_regex × use_cleri) combinations must agree
            base = read_frames(path("ct.xyz"); use_regex=true, use_cleri=true)
            for ur in (true, false), uc in (true, false)
                @test read_frames(path("ct.xyz"); use_regex=ur, use_cleri=uc) == base
            end
        end

        @testset "write format strings" begin
            frame = Dict{String,Any}(
                "N_atoms" => 1,
                "cell" => [5.44 0.0 0.0; 0.0 5.44 0.0; 0.0 0.0 5.44],
                "pbc" => [true, true, true],
                "info" => Dict{String,Any}(),
                "arrays" => Dict{String,Any}(
                    "species" => ["Si"],
                    "pos" => reshape([1/3, 2/3, 0.123456789123456789], 3, 1)))
            # 16 decimal places retains near-full double precision...
            write_frame(path("hp.xyz"), frame; fmt_f="%21.16f")
            hp = read_frame(path("hp.xyz"))["arrays"]["pos"]
            @test maximum(abs.(hp - frame["arrays"]["pos"])) <= 1e-15
            # ...while the default format (%16.8f) truncates at 8
            write_frame(path("lp.xyz"), frame)
            lp = read_frame(path("lp.xyz"))["arrays"]["pos"]
            @test 1e-12 < maximum(abs.(lp - frame["arrays"]["pos"])) <= 1e-8
        end
    end
end
