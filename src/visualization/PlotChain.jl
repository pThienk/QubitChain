# Collection of plotting interface fucntions and its component calls

# Includes the plot recipes
include("PlotRecipes.jl")

function plot_chain(chain_data::Union{QChain, QChainData}; filename::String="", settings...)
    qubits_plot_data::Vector{Tuple}, N::Int64, boundary_cond::Symbol = plot_obj_convert(chain_data, :graph)

    
end

function add_qubit_plot!(figure::Figure, plot_data::Vector{Tuple}, graph_num::Int64=10; sub_plot_settings...)
    
end

function plot_parameters()
    
end

function bloch_chain_anim()

end