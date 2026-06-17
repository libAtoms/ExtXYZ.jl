using ExtXYZ
using Test

# Error-message propagation and EOF detection (libextxyz >= 0.4 returns a
# diagnostic message; ExtXYZ.jl must distinguish end-of-file from parse errors)
@testset "Errors and EOF" begin
    good_frame(x) = """2
Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44" Properties=species:S:1:pos:R:3
Si 0.0 0.0 $x
Ge 1.36 1.36 1.36
"""
    mktempdir() do dir
        path(name) = joinpath(dir, name)

        @testset "parse errors raise descriptive exceptions" begin
            # non-integer natoms header: not whitespace, so NOT silent EOF
            write(path("badnat.xyz"), "abc\nfoo\n")
            err = @test_throws ErrorException read_frames(path("badnat.xyz"))
            @test occursin("Failed to parse int natoms", err.value.msg)

            # unterminated quoted string in the comment line
            write(path("badcomment.xyz"),
                  "2\nLattice=\"5.44 0.0 0.0 Properties=species:S:1:pos:R:3\nSi 0.0 0.0 0.0\nGe 1.0 1.0 1.0\n")
            err = @test_throws ErrorException read_frames(path("badcomment.xyz"))
            @test occursin("extxyz parse error", err.value.msg)

            # atom line with too few columns for the declared Properties
            write(path("badcols.xyz"),
                  "2\nLattice=\"5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44\" Properties=species:S:1:pos:R:3\nSi 0.0 0.0\nGe 1.0 1.0 1.0\n")
            err = @test_throws ErrorException read_frames(path("badcols.xyz"))
            @test occursin("atom line", err.value.msg)
            err = @test_throws ErrorException read_frames(path("badcols.xyz"); use_regex=false)
            @test occursin("atom line", err.value.msg)

            # read_frame goes through the same path
            @test_throws ErrorException read_frame(path("badcols.xyz"))
        end

        @testset "EOF detection" begin
            body = join(good_frame.(1:10))

            # file ending exactly at the final newline
            write(path("exact.xyz"), body)
            @test length(read_frames(path("exact.xyz"))) == 10

            # trailing blank lines
            write(path("blank.xyz"), body * "\n\n\n")
            @test length(read_frames(path("blank.xyz"))) == 10

            # trailing whitespace-only lines
            write(path("ws.xyz"), body * "   \n\t\n  ")
            @test length(read_frames(path("ws.xyz"))) == 10

            # iread_frames terminates cleanly at EOF
            n = 0
            for frame in iread_frames(path("exact.xyz"))
                n += 1
            end
            @test n == 10

            # a truncated final frame (header claims 5 atoms, only 3 present)
            # is dropped silently; preceding complete frames are returned
            write(path("truncated.xyz"), join(good_frame.(1:3)) * """5
Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44" Properties=species:S:1:pos:R:3
Si 0.0 0.0 0.0
Ge 1.0 1.0 1.0
Si 2.0 2.0 2.0
""")
            @test length(read_frames(path("truncated.xyz"))) == 3
        end
    end
end
