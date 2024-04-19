# Plot recipes for converting data structures to plottable arrays

function plot_obj_covert(chain_data::Union{QChain, QChainData})

    N::Int64 = chain_data.N
    boundary_cond::Symbol = chain_data.boundary_cond

    for i ∈ 1:N 
        
    end
    
    if chain_data isa QChain

        

    end
end

function convert_to_graph_tuple(qubit::QubitData)::Tuple
    return zip(qubit.t, abs.(qubit.dot_1) .^2, abs.(qubit.dot_2) .^2) |> collect
end