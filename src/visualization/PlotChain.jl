# Collection of plotting interface fucntions and its component calls

# Includes the plot recipes
include("PlotRecipes.jl")

function plot_chain(chain_data::Union{QChain, QChainData}; filename::String="", settings...)
    qubits_plot_data::Vector{Tuple}, N::Int64, boundary_cond::Symbol = plot_obj_convert(chain_data, :graph)
    acc_plot_num::Int64 = 0

    qubit_figs::Vector{Figure} = []
    param_fig::Figure = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (800, 700))

    add_param_plots!(param_fig, param_to_scatter_arr(chain_data), chain_data.boundary_cond)

    for i ∈ 1:floor(Int64, N/10) 
        push!(qubit_figs, Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size=(800, 1000)))
        add_qubit_plot!(qubit_figs[i], qubits_plot_data, acc_plot_num)
        acc_plot_num += 10
    end

    if N % 10 > 0
        push!(qubit_figs, Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size=(800, 1000)))
        add_qubit_plot!(qubit_figs[i], qubits_plot_data, acc_plot_num, N % 10)
    end
    
end

function add_qubit_plot!(figure::Figure, plot_data::Vector{Tuple}, plot_ind::Int64, graph_num::Int64=10)
    
end

function add_param_plots!(figure::Figure, param_plot_data::Tuple, boundary_cond::Symbol)
    J, Δ, ϵ, (J_avg, Δ_avg, ϵ_avg) = param_plot_data

    J_grid = figure[1, 1] = GridLayout()
    Δ_grid = figure[2, 1] = GridLayout()
    ϵ_grid = figure[3, 1] = GridLayout()

    ax_g1 = Axis(J_grid, xlabel = L"\textrm{J}~(J)", title="Fluctuation of the J, Δ, ϵ Parameters")
    ax_g2 = Axis(Δ_grid, xlabel = L"Δ~(J)")
    ax_g3 = Axis(ϵ_grid, xlabel = L"ϵ~(J)", ylabel = L"N")

    if boundary_cond == :open_ended
        scatter!(ax_g1, zip(J, 1:ndims(J) |> collect)[1:end-1] |> collect, label="J", color=:purple)
    else
        scatter!(ax_g1, zip(J, 1:ndims(J) |> collect) |> collect, label="J", color=:purple)
    end
    
    scatter!(ax_g2, zip(Δ, 1:ndims(Δ) |> collect) |> collect, label="Δ", color=:blue)
    scatter!(ax_g3, zip(ϵ, 1:ndims(ϵ) |> collect) |> collect, label="ϵ", color=:yellow)

    lines!(ax_g1, 0..20, x -> J_avg, linestyle=:dash, color=:black, alpha=0.75)
    lines!(ax_g2, 0..20, x -> Δ_avg, linestyle=:dash, color=:black, alpha=0.75)
    lines!(ax_g3, 0..20, x -> ϵ_avg, linestyle=:dash, color=:black, alpha=0.75)

    resize_to_layout!(figure)
end

function bloch_chain_anim()

end