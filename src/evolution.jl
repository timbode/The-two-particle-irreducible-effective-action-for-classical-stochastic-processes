# vertical update
function update!(i::Int, j::Int, solver::Function, M::Model)
    foreach(zip(M.G.data, M.rhs.vert(M))) do (g, f)
        for b in 1:M.P.dim, a in 1:M.P.dim
            g[a, b][i + 1, j] = g[a, b][i, j] + M.P.delta_t * solver(i, i -> f(a, b, i, j))
        end
    end        
end

# diagonal update
function update!(i::Int, solver::Function, M::Model)
    foreach(zip(M.G.data, M.rhs.diag(M))) do (g, f)
        for b in 1:M.P.dim, a in 1:M.P.dim
            g[a, b][i + 1, i + 1] = g[a, b][i, i] + M.P.delta_t * solver(i, i -> f(a, b, i, i))
        end
    end        
end

# --------------------------------------------------------------------------------------------------

# symmetry operations
function reflect!(i::Int, j::Int, M::Model)
    for g in M.G.data
        for b in 1:M.P.dim, a in 1:M.P.dim
            # classical Green functions
            g[a, b][j, i+1] = g[b, a][i+1, j]
        end
    end     
end

# --------------------------------------------------------------------------------------------------

# vertical memory
function memorize!(i::Int, j::Int, M::Model) 
    # compute right-hand sides
    foreach(zip(M.M.data, rhs_vert(M))) do (m, f)
        for b in 1:M.P.dim, a in 1:M.P.dim
            m[a, b][i, j] = f(a, b, i, j)
        end
    end     
end

# diagonal memory
function memorize!(i::Int, M::Model)
    # compute right-hand sides
    foreach(zip(M.M.diagonal, rhs_diag(M))) do (m, f)
        for b in 1:M.P.dim, a in 1:M.P.dim
            m[a, b][i] = f(a, b, i, i)
        end
    end     
end

# compute right-hand sides above diagonal where needed
# (which is for points up to three steps below the diagonal)
function memorize_above_diag!(i::Int, j::Int, M::Model)
    for k in 1:(3 - (i - j))
        memorize!(i - k, j, M)
    end
end

# --------------------------------------------------------------------------------------------------

function two_time_step(i_start::Int, i_end::Int, 
                       vert_predict::Function, 
                       diag_predict::Function,
                       vert_correct::Function,
                       diag_correct::Function)
    
    for i in i_start:i_end # t-loop
        
        map(j -> vert_predict(i, j), 1:i) # t'-loop
        
        # get next diagonal element
        diag_predict(i)
        
        map(j -> vert_correct(i, j), 1:i) # t'-loop     

        # get next diagonal element
        diag_correct(i)      
        
    end
end

function evolution!(M::Model)
    
    # take a predictor step to (i+1, j) and reflect
    vert_predict = (i, j) -> (update!(i, j, M.solver.predictor, M); 
                              reflect!(i, j, M)) 

    diag_predict =  i     ->  update!(i, M.solver.predictor, M)

    # take a corrector step to (i+1, j) and reflect
    vert_correct = (i, j) -> (update!(i, j, M.solver.corrector, M); 
                              reflect!(i, j, M))

    diag_correct =  i     ->  update!(i, M.solver.corrector, M)
    
    if M.P.initialize
        # Euler-Heun method
        two_time_step(1, 3, vert_predict, diag_predict, vert_correct, diag_correct)          
        
        # switch from Euler-Heun to ABM
        M = Model(M.P, M.G, M.M, Solver(AdamsBashforth, AdamsMoulton), M.rhs)
        
        # ABM method
        two_time_step(4, M.P.n, vert_predict, diag_predict, vert_correct, diag_correct)          
        
    else
        # ABM method with given memory
        two_time_step(4, M.P.n, 
            (i, j) -> (memorize!(i, j, M); 
                       memorize_above_diag!(i, j, M);
                       vert_predict(i, j)),      
            
                 i -> (memorize!(i, M);
                       diag_predict(i)),
            
            (i, j) -> (memorize!(i+1, j, M);
                       vert_correct(i, j)),
            
                 i -> (memorize!(i+1, M);
                       diag_correct(i))) 
    end
end