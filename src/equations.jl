# 2PI effective action for SDE with quadratic non-linearity and additive noise

# # 2PI-I approximation from the paper
# withIntegrals = false

# 2PI-2 approximation from the paper
withIntegrals = true

# ============================================================================================
# helper functions

# trapezoid integration
function ∫dt(k_a::Int, k_b::Int, F::Function)::ComplexF64
    return sum([F(k) for k in k_a:k_b]) - 0.5 * (F(k_a) + F(k_b))
end

# Heaviside
function Θ(t::Int64)::Int64
    return t > 0
end

# ============================================================================================
# rhs of equation for mean-field

function f_x(a::Int, b::Int, i::Int, j::Int, M::Model)::Complex{Float64}
    return 0.0im
end

function f_x_diag(a::Int, b::Int, i::Int, j::Int, M::Model)::Complex{Float64}
    x = M.G.data[1][1, 1]
    F = M.G.data[2][1, 1]
    return M.P.params[1] * x[i, i] + M.P.params[3] * (x[i, i]^2 + F[i, i])
end

# ============================================================================================
# rhs of two-time functions

# rhs of equation for statistical propagator
function f_F(a::Int, b::Int, i::Int, j::Int, M::Model)::Complex{Float64}
    x = M.G.data[1][1, 1]
    F = M.G.data[2][1, 1]
    K = M.G.data[3][1, 1]
    
    I = 0.0im
    if withIntegrals
        I = (M.P.delta_t * ∫dt(1, i, k -> 8M.P.params[3]^2 * K[i, k] * F[i, k] * F[k, j])
           + M.P.delta_t * ∫dt(1, j, k -> 4M.P.params[3]^2 * F[i, k]^2 * K[k, j]))
    end
    
    return (M.P.params[1] + 2 * M.P.params[3] * x[i, i]) * F[i, j] + M.P.params[2] * Θ(j - i) * K[i, j] + I/2
end

function f_F_diag(a::Int, b::Int, i::Int, j::Int, M::Model)::Complex{Float64}   
    x = M.G.data[1][1, 1]
    F = M.G.data[2][1, 1]    
    K = M.G.data[3][1, 1]  
    
    I = 0.0im
    if withIntegrals
        I = M.P.delta_t * ∫dt(1, i, k -> 24M.P.params[3]^2 * K[i, k] * F[i, k]^2)
    end    
    
    return (2 * M.P.params[1] + 4 * M.P.params[3] * x[i, i]) * F[i, i] + M.P.params[2] + I/2
end

# rhs of equation for response kernel
function f_K(a::Int, b::Int, i::Int, j::Int, M::Model)::Complex{Float64}
    x = M.G.data[1][1, 1]    
    F = M.G.data[2][1, 1]
    K = M.G.data[3][1, 1]    
    
    I = 0.0im
    if withIntegrals
        I = 4 * M.P.params[3]^2 * M.P.delta_t * ∫dt(i, j, k ->  F[i, k] * K[i, k] * K[k, j])
    end
    
    return (M.P.params[1] + 2 * M.P.params[3] * x[i, i]) * K[i, j] - I
end

function f_K_diag(a::Int, b::Int, i::Int, j::Int, M::Model)::Complex{Float64}   
    return 0.0im
end

# ============================================================================================
# package the rhs and define analogous functions
# accessing the already available memory

function rhs_vert(M::Model)::Array{Function, 1}
    return [(a, b, i, j) -> f_x(a, b, i, j, M), (a, b, i, j) -> f_F(a, b, i, j, M), (a, b, i, j) -> f_K(a, b, i, j, M)]
end

function rhs_diag(M::Model)::Array{Function, 1}
    return [(a, b, i, j) -> f_x_diag(a, b, i, j, M), (a, b, i, j) -> f_F_diag(a, b, i, j, M), (a, b, i, j) -> f_K_diag(a, b, i, j, M)]
end

function rhs_memory_vert(M::Model)::Array{Function, 1}
    return [(a, b, i, j) -> d[a, b][i, j] for d in M.M.data]
end

function rhs_memory_diag(M::Model)::Array{Function, 1}
    return [(a, b, i, j) -> d[a, b][i] for d in M.M.diagonal]
end
