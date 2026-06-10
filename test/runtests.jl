using Test

@testset "ExtXYZ.jl" begin
    include("dict.jl")
    include("errors.jl")
    include("strings.jl")
    include("stress.jl")
    include("atomsbase.jl")
end
