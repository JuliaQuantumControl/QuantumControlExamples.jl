using Test
using SafeTestsets
using Plots

unicodeplots()

# Note: comment outer @testset to stop after first @safetestset failure
@time @testset verbose = true "QuantumControlExample" begin

    @test true

    println("")

end
