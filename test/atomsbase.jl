using AtomsBase
using AtomsBaseTesting
using ExtXYZ
using Test
using Unitful
using UnitfulAtomic

@testset "Conversion AtomsBase -> Atoms" begin
    system = make_test_system().system
    test_approx_eq(system, Atoms(system))
end

@testset "Conversion AtomsBase -> dict (velocity)" begin
    system, atoms, atprop, sysprop, box, bcs = make_test_system()
    atoms = ExtXYZ.write_dict(Atoms(system))

    cell = zeros(3, 3)
    for i in 1:3
        cell[i, :] = ustrip.(u"Å", box[i])
    end

    @assert bcs == [Periodic(), Periodic(), DirichletZero()]
    @test atoms["pbc"]     == [true, true, false]
    @test atoms["N_atoms"] == 5
    @test atoms["cell"]    == cell

    info = atoms["info"]
    @test sort(collect(keys(info))) == ["charge", "extra_data", "multiplicity"]
    @test info["charge"]       == ustrip.(u"e_au", sysprop.charge)
    @test info["extra_data"]   == sysprop.extra_data
    @test info["multiplicity"] == sysprop.multiplicity

    arrays = atoms["arrays"]
    @test arrays["Z"]          == atprop.atomic_number
    @test arrays["species"]    == string.(atprop.atomic_symbol)
    @test arrays["mass"]       == ustrip.(u"u",  atprop.atomic_mass)
    @test arrays["pos"]        ≈  ustrip.(u"Å",  hcat(atprop.position...)) atol=1e-10
    @test arrays["velocities"] ≈  ustrip.(sqrt(u"eV"/u"u"),
                                          hcat(atprop.velocity...)) atol=1e-10

    expected_atkeys = ["Z", "charge", "covalent_radius", "magnetic_moment",
                       "mass", "pos", "species", "vdw_radius", "velocities"]
    @test sort(collect(keys(arrays))) == expected_atkeys
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
    @test atomic_number(atoms)    == [1]
    @test iszero(velocity(atoms, 1))
end

@testset "Warning about setting invalid data" begin
    system = make_test_system(; extra_sysprop=(md=3u"u", symboldata=:abc),
                                extra_atprop=(massdata=3ones(5)u"u",
                                              atomic_symbol=[:D, :H, :C, :N, :He])).system
    atoms = @test_logs((:warn, r"Unitful quantity massdata is not yet"),
                       (:warn, r"Writing quantities of type Symbol"),
                       (:warn, r"Unitful quantity md is not yet"),
                       (:warn, r"Mismatch between atomic numbers and atomic symbols"),
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
    test_approx_eq(system, io_system; rtol=1e-6)
end
