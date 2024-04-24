# Interfacing funtions for connecting user calls to backend

"""
    Hamiltonian with Random Fluctuation

    Takes in median J0 and Δ0 values and standard interaction model (nearest neigbors only).
    Returns a vector of QHamiltonian_S for the specified model.
    Optional Args: J Fluctuation (in percentage, default: 10%), Δ Fluctuation (in percentage, default: 10%), Type of random distribution (default: Normal)

"""
function H_rand_model(J0::Float64, Δ0::Float64, ϵ::Float64, N::Int64, boundary_cond::Symbol, interaction_model::BitMatrix; distribution::Symbol=:uniform,
    J_fluc::Float64=0.1, Δ_fluc::Float64=0.1)::Vector{QHamiltonian_S}

    # psudorandom generator
    rand_algor = MersenneTwister(rand_seed)

    param_rand_pair::Vector{Float64} = [J0, Δ0]
    edge_interaction::BitMatrix = [0 0; 0 0]
    chain_H::Vector{QHamiltonian_S} = []
    
    for i ∈ 1:N
        if distribution == :uniform
            param_rand_pair = (rand(rand_algor, 2) .- 0.5) .* 2.0
            if boundary_cond == :open_ended && (i == 1 || i == N)
                edge_interaction = (i==1) ? [0 0; 0 1] : [1 0; 0 0]
                push!(chain_H, QHamiltonian_S(J = J0*(1 + J_fluc*param_rand_pair[1]), Δ = Δ0*(1 + Δ_fluc*param_rand_pair[2]), ϵ = ϵ,
                interaction_graph = edge_interaction))
            else
                push!(chain_H, QHamiltonian_S(J = J0*(1 + J_fluc*param_rand_pair[1]), Δ = Δ0*(1 + Δ_fluc*param_rand_pair[2]), ϵ = ϵ,
                interaction_graph = interaction_model))
            end
        elseif distribution == :normal
            dist_J = Normal(J0, J_fluc*J0/3)
            dist_Δ = Normal(Δ0, Δ_fluc*Δ0/3)

            param_rand_pair = [rand(rand_algor, dist_J), rand(rand_algor, dist_Δ)]
            if boundary_cond == :open_ended && (i == 1 || i == N)
                edge_interaction = (i==1) ? [0 0; 0 1] : [1 0; 0 0]
                push!(chain_H, QHamiltonian_S(J = param_rand_pair[1], Δ = param_rand_pair[2], ϵ = ϵ,
                interaction_graph = edge_interaction))
            else
                push!(chain_H, QHamiltonian_S(J = param_rand_pair[1], Δ = param_rand_pair[2], ϵ = ϵ,
                interaction_graph = interaction_model))
            end
        end
    end
    
    return chain_H
end

"""
    Initial Qubit Chain Builder

    Takes an existing QChain object and an undertermined number of range models
    Returns a QChainData encoded with the initial conditions for each qubits in the chain
        
"""
function initial_chain_model(N::Int64, model_range::Tuple{UnitRange{Int64}, Tuple{ComplexF64, ComplexF64}}...)::QChainInitial
    
    initial_cond_array::Vector{Tuple{ComplexF64, ComplexF64}} = []

    for i ∈ 1:N 
        for model ∈ model_range
            if i ∈ model[1]
               push!(initial_cond_array, model[2]) 
            end
        end
    end

    return QChainInitial(N = N, initial_states = initial_cond_array)
end

# Add quailty of life for users
"""
    Simulation Call

    Simulate the provided QChain object and mutates the object to store data within
    Before simulation, it prints out a summary of the system suppiled and ask for confirmation
    Data is not stored and the operation is aborted safely, if the simulation fails

"""
function simulate!(q_chain::QChain, initial_cond::QChainInitial, t_inv::Tuple{Float64, Float64}; solver_settings::NamedTuple=(), ode_settings::NamedTuple=())
    chain_DEQ_prob!, param_array::Vector{Float64} = parse_chain(q_chain)

    @info "Built DEQ system successfully"

    initial_array::Vector{ComplexF64} = I_extract_initial_cond(initial_cond) 

    ODE_prob = ODEProblem(chain_DEQ_prob!, initial_array, t_inv, param_array; ode_settings...)

    @info """ DEQ problem created, ready to simulate chain! Proceed? (Type "yes" to Proceed) """
    if readline() ≠ "yes"
        @warn "Aborting simulation..."
        return
    end

    try
        @time solution = solve(ODE_prob; solver_settings...)
    catch err
        @error "Simulation Failed"
        show(err)
        return
    end
    
    q_chain.full_sol = solution
    @info "ODESolution object successfully stored in the QChain object"

    for i ∈ 1:2:2q_chain.N 
        push!(q_chain.qubits, QubitData(dot_1_raw = solution[i, :], dot_2_raw = solution[i+1, :], dot_1 = solution(t_inv(1):t_step:t_inv(2))[i, :],
        dot_2 = solution(t_inv(1):t_step:t_inv(2))[i+1, :], t_inv(1):t_step:t_inv(2)))
    end
    @info "Interpolated data for step size $t_step successfully stored in the QChain object"

    @info "Simulation Complete"
end

function save_chain(q_chain::QChain; save::String="", inc_raw_t::Bool=true)
    chain_col_headers::Vector{String} = []
    param_col_headers::Vector{String} = ["J", "delta", "epsilon"]
    filename_chain::String = ""
    filename_param::String = ""
    chain_data::Vector = []

    push!(chain_col_headers, "t")
    push!(chain_data, q_chain.qubits[1].t)
    
    for i ∈ 1:q_chain.N 
        push!(chain_col_headers, string(i) * "_dot_1", string(i) * "_dot_2")
        push!(chain_data, q_chain.qubits[i].dot_1, q_chain.qubits[i].dot_2)
    end

    for i ∈ 1:q_chain.N 
        push!(chain_col_headers, string(i) * "_dot_raw_1", string(i) * "_dot_raw_2")
        push!(chain_data, q_chain.qubits[i].dot_1_raw, q_chain.qubits[i].dot_2_raw)
    end

    if inc_raw_t
        push!(chain_col_headers, "t_raw", "boundary_cond")
        push!(chain_data, q_chain.full_sol.t, [string(q_chain.boundary_cond)])
    else
        push!(chain_col_headers, "boundary_cond")
        push!(chain_data, [string(q_chain.boundary_cond)])
    end

    if save != ""
        filename_chain = string(split(save)[1]) * "_qubits" * ".csv"
        filename_param = string(split(save)[1]) * "_param" * ".csv"
    else
        id::String = string(rand(1:1000000))
        filename_chain = "qubit_chain_sim_id-" * id * ".csv"
        filename_param = "qubit_chain_param_id-" * id * ".csv"
    end

    write_to_csv(filename_chain, chain_col_headers, chain_data...)
    write_to_csv(filename_param, param_col_headers, parse_chain(q_chain; op_type=:param_arr)...)
end

function load_data(filename::String, col_headers::String...; ret_type::Symbol=:qchain)::Union{QChainData, DataFrame}
    data_table::DataFrame = read_from_csv(filename)
    
    if isempty(col_headers)
        if ret_type == :dataframe
            return data_table
        elseif ret_type == :qchain
            return encode_data(data_table, ret_type)
        else
            return data_table
        end
    else
        if ret_type == :array
            return encode_data(data_table, ret_type, col_headers...)
        else
            return data_table
        end
    end
end

function load_parameters(filename::String, col_headers::String...; ret_type::Symbol=:array)::Union{Tuple, DataFrame}
    data_table::DataFrame = read_from_csv(filename)

    if !isempty(col_headers)
        if ret_type == :dataframe
            return data_table
        elseif ret_type == :array
            return encode_data(data_table, ret_type, col_headers...)
        end
    else
        return data_table
    end
end

function visualize(q_chain::Union{QChain, QChainData}; type::Symbol=:graph, save::String="", plot_settings...)

    if type == :graph
        plot_chain(q_chain; filename=string(split(save, ".")[1]), plot_settings...)
    elseif type == :anim

    else
        throw("Invalid plot option!")
    end

end

function set_solution_step_size(δt::Float64)
    global t_step = δt
end

function set_rand_seed(seed::Int64)
    global rand_seed = seed
end

function set_t_ticks(tick::Float64)
    global t_ticks = tick
end

function set_param_y_ticks(tick::Float64)
    global param_y_ticks = tick
end