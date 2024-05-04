# Core functionalities of the program

# Constants
const σx::Matrix{Complex{Float64}} = [0 1; 1 0]

const σy::Matrix{Complex{Float64}} = [0 -im; im 0]

const σz::Matrix{Complex{Float64}} = [1 0; 0 -1]

const h::Float64 = 6.62607015e-34 # J ∘ S

const ħ::Float64 = 1.0 # h/2π J ∘ S

# Indexing
const PARAM_RET_IND::Int64 = 1

const DIM_RET_IND::Int64 = 2

const INTER_GRAPHS_RET_IND::Int64 = 3

const J_RET_IND::Int64 = 4

const Δ_RET_IND::Int64 = 5

const ϵ_RET_IND::Int64 = 6

# Initial state types
const CAT::Tuple{ComplexF64, ComplexF64} = (1/√2, -1/√2)

const ANTI_CAT::Tuple{ComplexF64, ComplexF64} = (1/√2, 1/√2)

const Y_PLUS::Tuple{ComplexF64, ComplexF64} = (1/√2, 1/√2im)

const Y_MINUS::Tuple{ComplexF64, ComplexF64} = (1/√2, -1/√2im)

const DOT_0::Tuple{ComplexF64, ComplexF64} = (1, 0)

const DOT_1::Tuple{ComplexF64, ComplexF64} = (0, 1)

# Global variables
t_step::Float64 = 1e-3

t_ticks::Float64 = 0.1

param_y_ticks::Float64 = 0.05

rand_seed::Int64 = 1876

# Units
time_unit = u"10^-16 * s"


function parse_chain(q_chain::QChain; op_type::Symbol=:build_ode)::Tuple
    param_collection::Vector{Float64}, N::Int64, interaction_graphs::Vector{BitMatrix}, c_J::Vector{Float64}, c_Δ::Vector{Float64},
    c_ϵ::Vector{Float64} = I_extract_param_H(q_chain.H)

    if op_type == :build_ode
        ODE_prob! = build_DEQ_sys(N, q_chain.boundary_cond, interaction_graphs)

        return (ODE_prob!, param_collection)
    else # op_type == :param_arr
        return (c_J, c_Δ, c_ϵ) # return parameter arrays
    end
    
end

function I_extract_initial_cond(initial_cond::QChainInitial)::Vector{ComplexF64}
    initial_array::Vector{ComplexF64} = []

    for i ∈ eachindex(initial_cond.initial_states) 
        push!(initial_array, initial_cond.initial_states[i]...)
    end

    return initial_array
end

function I_extract_param_H(c_hamiltonian::Vector{QHamiltonian_S})::Tuple
    c_J::Vector{Float64} = []
    c_Δ::Vector{Float64} = []
    c_ϵ::Vector{Float64} = []

    param_collection::Vector{Float64} = []
    interaction_graphs::Vector{BitMatrix} = []
    for i ∈ eachindex(c_hamiltonian) 
        push!(c_J, c_hamiltonian[i].J)
        push!(c_Δ, c_hamiltonian[i].Δ)
        push!(c_ϵ, c_hamiltonian[i].ϵ)
        push!(param_collection, c_hamiltonian[i].J, c_hamiltonian[i].Δ, c_hamiltonian[i].ϵ)
        push!(interaction_graphs, c_hamiltonian[i].interaction_graph)
    end

    return (param_collection, length(c_hamiltonian), interaction_graphs, c_J, c_Δ, c_ϵ)
end

function encode_data(data::DataFrame, type::Symbol, col_headers::String...)::Union{QChainData, Tuple}
    
    if type == :qchain
        t::Vector{Float64} = data[1:end, :t]
        t_raw::Vector{Float64} = (data[2:end, :t_raw] |> skipmissing |> collect)
        qubits::Vector{QubitData} = []
        N::Int64 = (length(data[1, :]) - 3) ÷ 4
        if N % 2 == 0
            
        end
        for i ∈ eachindex(data[1, :])[1:2:2N-2]
            push!(qubits, QubitData(dot_1 = parse.(ComplexF64, data[1:end, i+1]), dot_2 = parse.(ComplexF64, data[1:end, i+2]),
            dot_1_raw = parse.(ComplexF64, (data[1:end, 2N + 1 + i] |> skipmissing |> collect)),
            dot_2_raw = parse.(ComplexF64, (data[1:end, 2N + 2 + i] |> skipmissing |> collect), t = t)))
        end

        return QChainData(qubits = qubits, N = N, t_raw = t_raw, boundary_cond = Symbol(data[1, end]))
    elseif type == :array
        if isempty(col_headers)
            throw("Call must actually contains columns!")
        end
        extracted_data::Vector = []
         
        for col ∈ Symbol.(col_headers)
            push!(extracted_data, (data[1:end, col] |> skipmissing |> collect))
        end

        return Tuple(extracted_data)
    end
end

function ⊗(A::Matrix{String}, B::Vector{String})::Vector{String}
    result_vec::Vector{String} = []
    for i ∈ eachindex(B[:])
        column_str::Vector{String} = []
        for j ∈ eachindex(A[1, :]) 
            if A[i, j] ≠ "" && B[j] ≠ ""
                if A[i, j] == "1.0 + 0.0im"
                    push!(column_str, B[j])
                elseif B[j] == "1.0 + 0.0im"
                    push!(column_str, A[i, j])
                elseif B[j] == "-1.0 + 0.0im"
                    push!(column_str, join(["-", A[i, j]]))
                elseif A[i, j] == "-1.0 + 0.0im"
                    push!(column_str, join("-", B[j]))
                else
                    push!(column_str, join(["(", A[i, j], ")*(", B[j], ")",]))
                end
            end
        end
        if column_str == []
            push!(result_vec, "")
        else
            push!(result_vec, join(column_str, " + "))
        end 
    end

    return result_vec
end

function ⊗(c::Union{Number, String}, A::Matrix{String})
    result_mat::Matrix{String} = ["" ""; "" ""]
    for i ∈ eachindex(A) 
        if A[i] == ""
            result_mat[i] = ""
        elseif A[i] == "1.0 + 0.0im"
            result_mat[i] = join([c,""])
        elseif A[i] == "-1.0 + 0.0im"
            result_mat[i] = join(["-", c,""])
        else
            result_mat[i] = join(["(", c, ")", "*", "(", A[i], ")"])
        end
    end

    return result_mat
end

# Matrix-Matrix expression multiplication (For the future add-on)
function ⊗(A::Matrix{String}, B::Matrix{String})::Matrix{String}

    if length(A[1, :]) ≠ length(B[:, 1])
        throw(DimensionMismatch("Invalid matrix multiplication call"))
    end

    result_mat::Matrix{String} = parse_str_matrix(zeros(length(A[:, 1]), length(B[1, :])))
    
    for ind_A_1 ∈ A[:, 1]
        for ind_B_2 ∈ B[1, :] 
            for n ∈ B[:, 1] 
                if A[ind_A_1, n] ≠ "" && B[n, ind_B_2] ≠ ""
                    if A[i, j] == "1.0 + 0.0im"
                            push!(element_str, B[n, ind_B_2])
                    elseif B[n, ind_B_2] == "1.0 + 0.0im"
                            push!(element_str, A[ind_A_1, n])
                    elseif B[n, ind_B_2] == "-1.0 + 0.0im"
                            push!(element_str, join(["-", A[ind_A_1, n]]))
                    elseif A[ind_A_1, n] == "-1.0 + 0.0im"
                            push!(element_str, join("-", B[n, ind_B_2]))
                    else
                            push!(element_str, join(["(",  A[ind_A_1, n], ")*(", B[n, ind_B_2], ")",]))
                    end
                end
            end

            result_mat[ind_A_1, ind_B_2] = join(element_str, " + ")
        end
    end

    return result_mat
end

function ⊕(A::AbstractArray{String}, B::AbstractArray{String})::AbstractArray{String}
    sum_array::AbstractArray{String} = []
    if ndims(A) == ndims(B)
        for i in eachindex(A)
            if A[i] == ""
                push!(sum_array, B[i])
            elseif B[i] == ""
                push!(sum_array, A[i])
            else
                push!(sum_array, join([A[i], " + ", B[i]]))
            end
        end
        return reshape(sum_array, ndims(A), ndims(A))
    end

    throw(DimensionMismatch("Array A and B must have the same shape."))
end

function mask_interaction_param_S(q_num::Int64, interaction_graph::BitMatrix, max_q_num::Int64)::Matrix{String}
    interaction_array::Matrix{String} = [join(["J", choose_max_or_ind(max_q_num, q_num, :upper), "_", q_num, "*", "(abs(", "A", 
    choose_max_or_ind(max_q_num, q_num, :upper), "_1)^2)"])  join(["J", choose_max_or_ind(max_q_num, q_num, :upper),
     "_", q_num, "*", "(abs(", "A", choose_max_or_ind(max_q_num, q_num, :upper), "_0)^2)"]); join(["J", q_num, "_", choose_max_or_ind(max_q_num, q_num, :lower),
      "*", "(abs(", "A", choose_max_or_ind(max_q_num, q_num, :lower), "_1)^2)"])  join(["J", q_num, "_", choose_max_or_ind(max_q_num, q_num, :lower),
       "*", "(abs(", "A", choose_max_or_ind(max_q_num, q_num, :lower), "_0)^2)"])]

    for i ∈ eachindex(interaction_graph)
        if interaction_graph[i] == 0
            interaction_array[i] = ""
        end
    end
    
    return interaction_array
end

function parse_str_matrix(M::Matrix{<:Number})
    result_mat::Matrix{String} = ["" ""; "" ""]
    for i ∈ eachindex(M) 
        if M[i] == 0
            result_mat[i] = ""
        else
            result_mat[i] = join([M[i], ""])
        end
    end

    return result_mat
end

function choose_max_or_ind(N::Int64, q_num::Int64, option::Symbol)
    if q_num-1 == 0 && option == :upper 
        return N
    elseif q_num == N && option == :lower
        return 1
    elseif option == :upper
        return q_num-1
    else
        return q_num+1
    end
end