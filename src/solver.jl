# Four-step Adams-Bashforth-Moulton (ABM) predictor-corrector

# predictor
function AdamsBashforth(i::Int, f::Function)::Complex{Float64}
    return (1.0/24.0) * (55.0 * f(i) - 59.0 * f(i-1) + 37.0 * f(i-2) - 9.0 * f(i-3))
end

# corrector
function AdamsMoulton(i::Int, f::Function)::Complex{Float64}
    return (1.0/24.0) * (9.0 * f(i+1) + 19.0 * f(i) - 5.0 * f(i-1) + f(i-2))
end

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

# Heun method for starting off

# predictor
function Euler(i::Int, f::Function)::Complex{Float64}
    return f(i)
end

# corrector
function Heun(i::Int, f::Function)::Complex{Float64}
    return 0.5 * (f(i) + f(i+1))
end