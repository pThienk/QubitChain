# Plot recipes for converting data structures to plottable arrays

function plot_obj_convert(chain_data::Union{QChain, QChainData}, type::Symbol)::Tuple

    N::Int64 = chain_data.N
    boundary_cond::Symbol = chain_data.boundary_cond

    if type == :graph
        qubits_plot_data::Vector{Tuple} = []
        for i ∈ 1:N 
            push!(qubits_plot_data, convert_to_graph_tuple(chain_data.qubits[i]))
        end

        return (qubits_plot_data, N, boundary_cond)
    elseif type == :bloch
        qubits_bloch_data::Vector{Tuple} = []
        for i ∈ 1:N 
            push!(qubits_bloch_data, convert_to_bloch_points(chain_data.qubits[i]))
        end
        
        return qubits_bloch_data
    else
        throw("Invalid conversion type: $(string(type))")
    end


end

function convert_to_graph_tuple(qubit::QubitData)::Tuple
    return (qubit.t, abs.(qubit.dot_1) .^2, abs.(qubit.dot_2) .^2) 
end

"""
    Points on the surface of a unit sphere: ̂r = cos(ϕ)sin(θ)̂x + sin(ϕ)sin(θ)̂y + cos(θ)̂z ,
    where dot_1 = cos(θ) + sin(θ)im, dot_2 = cos(ϕ) + sin(ϕ)im
"""
function convert_to_bloch_points(qubit::QubitData)::Tuple
    z::Vector{Float64} = real.(qubit.dot_1)
    y::Vector{Float64} = imag.(qubit.dot_2) .* imag.(qubit.dot_1)
    x::Vector{Float64} = real.(qubit.dot_2) .* imag.(qubit.dot_1)

    return (x, y, z)
end

function param_to_scatter_arr(param_data::Union{QChain, Vector{QHamiltonian_S}})::Tuple
    
    if param_data isa QChain
        J, Δ, ϵ = I_extract_param_H(param_data.H)[J_RET_IND:ϵ_RET_IND]
        
        return (J, Δ, ϵ, param_data.N, (mean(J), mean(Δ), mean(ϵ)))
    else
        J, Δ, ϵ = I_extract_param_H(param_data)[J_RET_IND:ϵ_RET_IND]

        return (J, Δ, ϵ, ndims(param_data), (mean(J), mean(Δ), mean(ϵ)))
    end
end
