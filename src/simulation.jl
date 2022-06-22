function simulation(n::Int, T::Float64, dim::Int, params::Array{ComplexF64, 1}, initial_data::Array{Array{ComplexF64, 2}, 1})::Array{Array{Array{ComplexF64, 2}, 2}, 1}
    
    delta_t = T/n
    
    # fine-grain the initial Euler-Heun points 
    # to get better starting values for the multi-step method
    n_ini = Int(3 * 2^5)
    delta_t_ini = (3 * delta_t)/n_ini
    
    # initial Green functions
    data_ini = [[zeros(ComplexF64, n_ini + 1, n_ini + 1) for _ in 1:dim, _ in 1:dim] for _ in initial_data]

    # initial conditions
    for (k, d) in enumerate(initial_data)
        for b in 1:dim, a in 1:dim
            data_ini[k][a, b][1, 1] = d[a, b]
        end       
    end

    pars_ini = Parameters(true, n_ini, delta_t_ini, dim, params)
    G_ini = GreenFunctions(data_ini)
    memory_ini = Memory([[zeros(ComplexF64, 1, 1) for _ in 1:1, _ in 1:1], 
                         [zeros(ComplexF64, 1, 1) for _ in 1:1, _ in 1:1]], 
                        [[zeros(ComplexF64, 1) for _ in 1:1, _ in 1:1]])
    solver_ini = Solver(Euler, Heun)
    rhs_ini = RHS(rhs_vert, rhs_diag)    
    model_ini = Model(pars_ini, G_ini, memory_ini, solver_ini, rhs_ini)     
    evolution!(model_ini)
    
#     ================================================================================
    
    # Green functions
    data = [[zeros(ComplexF64, n + 1, n + 1) for _ in 1:dim, _ in 1:dim] for _ in initial_data]

    # memory
    memory = [[zeros(ComplexF64, n + 1, n + 1) for _ in 1:dim, _ in 1:dim] for _ in initial_data]
    memory_diag = [[zeros(ComplexF64, n + 1) for _ in 1:dim, _ in 1:dim] for _ in initial_data]

    # Setting the first four points to start-off the multi-step method
    # using the fine-grained initial data
    idx = Int(n_ini/3)
    for (k, d) in enumerate(data_ini)
        for j in 1:4, i in 1:4, b in 1:dim, a in 1:dim
            data[k][a, b][i, j] = d[a, b][(i - 1) * idx + 1, (j - 1) * idx + 1]
        end        
    end
    
    pars = Parameters(false, n, delta_t, dim, params)
    G = GreenFunctions(data)
    M = Memory(memory, memory_diag)
    solver = Solver(AdamsBashforth, AdamsMoulton)
    rhs = RHS(rhs_memory_vert, rhs_memory_diag)
    model = Model(pars, G, M, solver, rhs)

    # Setting the first four cache values for the multi-step method
    for j in 1:4, i in 1:4
        memorize!(i, j, model)
        memorize!(i, model)
    end    
    
    evolution!(model) 
  
    return data
end