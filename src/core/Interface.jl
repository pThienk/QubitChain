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
function initial_chain_model(N::Int64, model_range::Tuple{AbstractRange{Int64}, Tuple{ComplexF64, ComplexF64}}...)::QChainInitial
    
    # Possible modification to include support for Vector{ComplexF64} or Vector{Tuple{ComplexF64, ComplexF64}} 

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
function simulate!(q_chain::QChain, initial_cond::QChainInitial, t_inv::Tuple{Float64, Float64}; solver_settings=(), ode_settings=(), comf_bypass=false)
    chain_DEQ_prob!, param_array::Vector{Float64} = parse_chain(q_chain)

    @info "Built DEQ system successfully"

    initial_array::Vector{ComplexF64} = I_extract_initial_cond(initial_cond) 

    tspan::Tuple{Float64, Float64} = t_inv
    if t_inv[1] > 0.0
        tspan = (0.0, t_inv[2])
    end

    ODE_prob = ODEProblem(chain_DEQ_prob!, initial_array, tspan, param_array; ode_settings...)

    @info """ DEQ problem created, ready to simulate chain! Proceed? (Type "yes" to proceed) """
    if !comf_bypass && readline() ≠ "yes"
        @warn "Aborting simulation..."
        return
    end
    
    @time solution = solve(ODE_prob; solver_settings...)
    @info "Simulated successfully!"
    
    q_chain.full_sol = solution
    @info "ODESolution object successfully stored in the QChain object"

    for i ∈ 1:2:2q_chain.N 
        push!(q_chain.qubits, QubitData(dot_1_raw = solution[i, :], dot_2_raw = solution[i+1, :], dot_1 = solution(t_inv[1]:t_step:t_inv[2])[i, :],
        dot_2 = solution(t_inv[1]:t_step:t_inv[2])[i+1, :], t = t_inv[1]:t_step:t_inv[2]))
    end
    @info "Interpolated data for step size $t_step successfully stored in the QChain object"

    @info "Simulation Complete"
end

"""
    Save Chain

    Saves data within a QChain object to a .csv file
    Args: QChain q_chain
    Optional Keyword Args: String save (if not provided, the function will use default names with randomly generated ids)
                   Bool inc_raw_t (determines if the raw time from ODE Solver is included, default=true)

"""
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
        filename_chain = string(split(save, ".")[1]) * "_qubits" * ".csv"
        filename_param = string(split(save, ".")[1]) * "_param" * ".csv"
    else
        id::String = string(rand(1:1000000))
        filename_chain = "qubit_chain_sim_id-" * id * ".csv"
        filename_param = "qubit_chain_param_id-" * id * ".csv"
    end

    write_to_csv(filename_chain, chain_col_headers, chain_data...)
    write_to_csv(filename_param, param_col_headers, parse_chain(q_chain; op_type=:param_arr)...)
    @info "Qubit dynamic data and parameters saved successfully!"
end

"""
    Load Qubit Data

    Loads qubit data from .csv file
    Args: String filename
    Optional Args: String... col_headers (let you choose specific columns by their names)
    Optional Keyword Args: Symbol ret_type (two options are :qchain and :dataframe, default=:qchain)

"""
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

"""
    Load Parameters Data

    Loads parameters data from .csv file
    Args: String filename
    Optional Args: String... col_headers (let you choose specific columns by their names)
    Optional Keyword Args: Symbol ret_type (two options are :array and :dataframe, default=:array)

"""
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

"""
    Visualize

    Visualizes either a QChain object or a QChainData object
    Args: Union{QChain, QChainData} q_chain
    Optional Keyword Args: Symbol type (type of visualization, two options are :graph and :anim (Not yet supported), default=:graph)
                           String save (filename for saving visual data, not saved if left empty)
                           plot_settings... (optional settings for the plotting functions)
    
"""
function visualize(q_chain::Union{QChain, QChainData}; type::Symbol=:graph, save::String="", param_plot_settings=(), qubit_plot_settings=())

    if type == :graph
        CairoMakie.activate!()
        plot_chain(q_chain; filename=string(split(save, ".")[1]), param_plot_settings=param_plot_settings, qubit_plot_settings=qubit_plot_settings)
    elseif type == :anim
        GLMakie.activate!()
        
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

function set_t_unit(unit::Quantity)
    global time_unit = unit
end