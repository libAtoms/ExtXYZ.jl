using AtomsBase
using AtomsBaseTesting
using ExtXYZ
using Test
using Unitful
using UnitfulAtomic
using AtomsBase: AbstractSystem


function simple_test_approx_eq(sys1, sys2; test_cell = true)
    @test all(position(sys1, :) .≈ position(sys2, :))
    @test all(velocity(sys1, :) .≈ velocity(sys2, :))
    @test mass(sys1, :) ≈ mass(sys2, :)
    @test species(sys1, :) == species(sys2, :)
    @test atomic_number(sys1, :) == atomic_number(sys2, :)
    @test atomic_symbol(sys1, :) == atomic_symbol(sys2, :)
    @test periodicity(sys1) == periodicity(sys2)
    @test all(cell_vectors(sys1) .≈ cell_vectors(sys2))
    if test_cell 
        @test cell(sys1) == cell(sys2)
    end
end    

@testset "Conversion AtomsBase -> Atoms" begin
    system = make_test_system().system
    simple_test_approx_eq(system, Atoms(system))
end

@testset "Conversion AtomsBase -> dict (velocity)" begin
    tmp = make_test_system()
    system = tmp.system
    atoms  = tmp.atoms
    atprop = tmp.atprop
    sysprop = tmp.sysprop
    box = tmp.cell.cell_vectors
    bcs = tmp.cell.periodicity
    atoms = ExtXYZ.write_dict(Atoms(system))

    c3ll = zeros(3, 3)
    for i in 1:3
        c3ll[i, :] = ustrip.(u"Å", box[i])
    end

    @assert bcs == (true, true, false)
    @test atoms["pbc"]     == [true, true, false]
    @test atoms["N_atoms"] == 5
    @test atoms["cell"]    == c3ll

    info = atoms["info"]
    @test sort(collect(keys(info))) == ["charge", "extra_data", "multiplicity"]
    @test info["charge"]       == ustrip.(u"e_au", sysprop.charge)
    @test info["extra_data"]   == sysprop.extra_data
    @test info["multiplicity"] == sysprop.multiplicity

    arrays = atoms["arrays"]
    @test arrays["Z"]          == AtomsBase.atomic_number.(atprop.species)
    @test arrays["species"]    == string.(atprop.species)
    # mass is not written to the dict because the mass == mass(element)
    # @test arrays["mass"]       == ustrip.(u"u",  atprop.atomic_mass)
    @test arrays["pos"]        ≈  ustrip.(u"Å",  hcat(atprop.position...)) atol=1e-10
    @test arrays["velocities"] ≈  ustrip.(sqrt(u"eV"/u"u"),
                                          hcat(atprop.velocity...)) atol=1e-10

    expected_atkeys = ["Z", "species", "charge", "covalent_radius", 
                       "magnetic_moment", "mass", "pos", "vdw_radius", "velocities"]
    @test sort(collect(keys(arrays))) == sort(expected_atkeys)
    @test arrays["magnetic_moment"] == atprop.magnetic_moment
    @test arrays["vdw_radius"]      == ustrip.(u"Å", atprop.vdw_radius)
    @test arrays["covalent_radius"] == ustrip.(u"Å", atprop.covalent_radius)
    @test arrays["charge"]          == ustrip.(u"e_au", atprop.charge)
end

@testset "Conversion AtomsBase -> dict (no velocity)" begin
    system = make_test_system(; drop_atprop=[:velocity]).system
    atoms  = ExtXYZ.write_dict(Atoms(system))
    @test iszero(atoms["arrays"]["velocities"])
end

@testset "Defaults dict -> Atoms" begin
    data = Dict("N_atoms" => 1, "pbc" => [true, true, true],
                "arrays" => Dict("pos" => [zeros(3)], "Z" => [1]),
                "cell" => randn(3, 3), "info" => Dict{String,Any}())
    atoms = Atoms(data)
    @test length(atoms) == 1
    @test all(periodicity(atoms))
    @test atomic_symbol(atoms, 1) == :H
    @test atomic_number(atoms, :)    == [1,]
    @test iszero(velocity(atoms, 1))
end

@testset "Warning about setting invalid data" begin
    system = make_test_system(; extra_sysprop=(md=3u"u", symboldata=:abc),
                                extra_atprop=(massdata=3ones(5)u"u",
                                              atomic_symbol=[:D, :H, :C, :N, :He])).system
    atoms = @test_logs((:warn, r"Unitful quantity massdata is not yet"),
                       (:warn, r"Writing quantities of type Symbol"),
                       (:warn, r"Unitful quantity md is not yet"),
                       match_mode=:any, ExtXYZ.write_dict(Atoms(system)))

    @test atoms["arrays"]["species"] == ["H", "H", "C", "N", "He"]
    @test atoms["arrays"]["Z"]       == [1, 1, 6, 7, 2]
end

@testset "AtomsBase -> IO -> Atoms" begin
    system = make_test_system().system
    io_system = mktempdir() do dir
        outfile = joinpath(dir, "data.extxyz")
        ExtXYZ.save(outfile, system)
        ExtXYZ.load(outfile)::AbstractSystem
    end
    # test_approx_eq(system, io_system; rtol=1e-4)
    simple_test_approx_eq(system, io_system)
end

@testset "Extra variables for atoms" begin
    text = """2
Lattice="2.614036117884091 0.0 0.0 0.0 2.6528336296738044 0.0 0.0 0.0 3.8250280122051756" Properties=species:S:1:pos:R:3:force:R:3:tags:I:1 config_type=FLD_TiAl spacegroup="P 1" virial="5.072173561696366 0.1220123768779895 0.6518229755809941 0.1220123768779895 4.667636799854875 0.5969893898844183 0.6518229755809941 0.5969893898844183 4.700422750506493" energy=-1703.64063822 unit_cell=conventional n_minim_iter=2 pbc="T T T"
Ti       1.30924260       1.32316179       1.62637131       0.86219000       0.78737000       2.65969000        0
Al       0.11095015       0.09471147      -0.05013464      -0.86219000      -0.78737000      -2.65969000        1
"""
    infile = tempname()
    outfile = tempname()
    open(infile, "w") do io
        print(io, text)
    end
    frame = read_frame(infile)
    atoms = Atoms(frame)
    @test all( atoms[1][:force] .== [0.86219000, 0.78737000, 2.65969000] )
    @test all( atoms[2][:force] .== [-0.86219000, -0.78737000, -2.65969000] )
    @test atoms[1][:tags] == 0
    @test atoms[2][:tags] == 1
    ExtXYZ.save(outfile, atoms)
    new_atoms = ExtXYZ.load(outfile)
    @test all( atoms[1][:force] .== new_atoms[1][:force] )
    @test all( atoms[2][:force] .== new_atoms[2][:force] )
    @test atoms[1][:tags] == new_atoms[1][:tags]
    @test atoms[2][:tags] == new_atoms[2][:tags]
end

@testset "AtomsBase missing cell" begin
    fname = tempname()
    try
        text = """13
        Properties=species:S:1:pos:R:3:forces:R:3 energy=-66.79083251953125 pbc="T T T"
        C       -5.13553286       0.00000000       0.00000000       0.00006186      -0.00000000       0.00000000
        H       -5.72485781      -0.91726500       0.00000000      -0.00001280       0.00000756       0.00000000
        H       -5.72485781       0.91726500       0.00000000      -0.00001280      -0.00000756      -0.00000000
        C       -3.82464290       0.00000000       0.00000000      -0.00004186       0.00000000       0.00000000
        C       -2.55226111       0.00000000       0.00000000       0.00003090      -0.00000000      -0.00000000
        C       -1.27494299       0.00000000       0.00000000      -0.00009426       0.00000000      -0.00000000
        C        0.00000000       0.00000000       0.00000000       0.00000000      -0.00000000      -0.00000000
        C        1.27494299       0.00000000       0.00000000       0.00009426       0.00000000       0.00000000
        C        2.55226111       0.00000000       0.00000000      -0.00003090      -0.00000000       0.00000000
        C        3.82464290       0.00000000       0.00000000       0.00004186       0.00000000      -0.00000000
        H        5.72485781       0.00000000       0.91726500       0.00001280      -0.00000000      -0.00000756
        H        5.72485781       0.00000000      -0.91726500       0.00001280       0.00000000       0.00000756
        C        5.13553286       0.00000000       0.00000000      -0.00006186      -0.00000000      -0.00000000
        """
        open(fname, "w") do io
            print(io, text)
        end

        atoms = ExtXYZ.load(fname)
        box = cell_vectors(atoms)
        @test [box[dim][dim] for dim=1:3] == [Inf * u"bohr" for dim=1:3]
    finally
        isfile(fname) && rm(fname)
    end
end

# TODO: This test is failing because ExtXYZ doesn't store a general cell, but 
#       always stores the cell vectors. Because of this, the equality test 
#       cannot pass.
@testset "AtomsBase isolated system" begin
    hydrogen = isolated_system([
            :H => [0, 0, 0.]u"Å",
            :H => [0, 0, 1.]u"Å"
    ])

    fname = tempname()
    try
        ExtXYZ.save(fname, hydrogen)
        new_sys = ExtXYZ.load(fname)
        # test_approx_eq(hydrogen, new_sys; rtol=1e-4)
        # note that test_cell = false only removes the equality test for 
        # the cell object, it still tests equality of the cell vectors and pbc
        simple_test_approx_eq(hydrogen, new_sys; test_cell=false)
    finally
        isfile(fname) && rm(fname)
    end
end