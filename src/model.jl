# summary of the parameters
struct Parameters   
    initialize::Bool # whether to use Euler-Heun to start off or not
    n::Int64 # number of time steps
    delta_t::Float64 # step size
    
    dim::Int64 # number of modes
    params::Vector{ComplexF64} # model parameters
end

# structure holding the Green functions
struct GreenFunctions    
    data::Vector{Matrix{Matrix{ComplexF64}}}
end

# structure holding the already evaluated right-hand sides
struct Memory    
    data::Vector{Matrix{Matrix{ComplexF64}}}
    diagonal::Vector{Matrix{Vector{ComplexF64}}}
end

# the solvers currently in use
# either Euler/Heun OR AdamsBashforth/AdamsMoulton
struct Solver
    predictor::Function
    corrector::Function
end

# the right-hand sides currently in use
# either the actual equations OR the already memorized right-hand sides
struct RHS
    vert::Function
    diag::Function
end

# putting it all into a single structure
struct Model
    P::Parameters
    G::GreenFunctions
    M::Memory
    solver::Solver
    rhs::RHS
end