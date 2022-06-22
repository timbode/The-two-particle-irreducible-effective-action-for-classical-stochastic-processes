using DrWatson
@quickactivate "The Two-Particle Irreducible Effective Action for Classical Stochastic Processes"

include("../src/ClassicalFixedStepTwoTimeSolver.jl")
using Test, LinearAlgebra

RTOL = 1e-8

@testset begin
    # time parameters
    delta_t = 10.0/2^10

    n = 1128
    T = n * delta_t    

    # parameters
    α = 1.0
    D = 1.5
    β = 0.15
    params = [α, D + 0.0im, β]

    # initial conditions
    x₀ = -10 * ones(ComplexF64, 1, 1)
    F₀ = 2 * ones(ComplexF64, 1, 1)
    K₀ = ones(ComplexF64, 1, 1)  
    
    x, F, K = ClassicalFixedStepTwoTimeSolver.simulation(n, T, 1, params, [x₀, F₀, K₀])

    @test (x[1, 1] |> diag)[1:100:end] |> real ≈ [-10.0, -7.560216590747483, -6.919839725068339, -6.693662291396823, -6.605157993364855, -6.56892231150991, -6.553754409829172, -6.547333089203879, -6.544598699543093, -6.5434308022701435, -6.542931207190734, -6.542717326809103] rtol = RTOL
    @test (F[1, 1] |> diag)[1:100:end] |> real ≈ [2.0, 0.6206459629492002, 0.6754802859897912, 0.7414641607162767, 0.778546981695146, 0.7967982112168673, 0.805260310254078, 0.8090538422055272, 0.8107211794856515, 0.8114456511535051, 0.8117584042624961, 0.8118929376112644] rtol = RTOL
end