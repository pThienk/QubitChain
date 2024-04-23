# Collection of plotting interface fucntions and its component calls

# Includes the plot recipes
include("PlotRecipes.jl")

function plot_chain(chain_data::Union{QChain, QChainData}; filename::String="", settings...)
    qubits_plot_data::Vector{Tuple}, N::Int64, boundary_cond::Symbol = plot_obj_convert(chain_data, :graph)
    acc_plot_num::Int64 = 0

    set_theme(theme_latexfonts())

    qubit_figs::Vector{Figure} = []
    param_fig::Figure = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (800, 700))

    add_param_plots!(param_fig, param_to_scatter_arr(chain_data), boundary_cond)

    for i ∈ 1:floor(Int64, N/10) 
        push!(qubit_figs, Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size=(800, 1000)))
        add_qubit_plot!(qubit_figs[i], qubits_plot_data, acc_plot_num)
        acc_plot_num += 10
    end

    if N % 10 > 0
        push!(qubit_figs, Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size=(800, 1000)))
        add_qubit_plot!(qubit_figs[i], qubits_plot_data, acc_plot_num, N % 10)
    end

    @info "Finished creating figures. \nDisplaying..."

    display_fig(param_fig)
    display_fig(qubit_figs)
    
    if filename ≠ ""
        save(filename * "_param.jpg", param_fig)

        for i ∈ eachindex(qubit_figs)
            save(filename * "_figure_" * string(('a' + i - 1)), qubit_figs[i])
        end

        @info "Saving completed."
    end
end

function add_qubit_plot!(figure::Figure, plot_data::Vector{Tuple}, plot_ind::Int64, graph_num::Int64=10)
    
    for i ∈ 1:graph_num

        t::Vector{Float64}, dot_1_data::Vector{Float64}, dot_2_data::Vector{Float64} = plot_data[i + plot_ind]

        grid_lvl = figure[i, 1:2] = GridLayout()
        grid_left = grid_lvl[1, 1] = GridLayout()
        grid_right = grid_lvl[1, 2] = GridLayout()

        ax_left = Axis(grid_left, xlabel = L"t~(ms)", ylabel = latexstring("qubit~$(plot_ind + i)"))
        ax_right = Axis(grid_right)
        lines!(ax_left, t, dot_1_data, label="dot 0")
        lines!(ax_left, t, dot_2_data, label="dot 1")

        if i ≠ graph_num
            hidexdecorations!(ax_left)
        end

        Legend(grid_right, ax_right)
        
        colgap!(grid_lvl, 10)
        rowgap!(grid_lvl, 5)
    end

    Label(figure[Top()], latexstring("Dynamic~of~qubit~$(plot_ind + 1)-$(plot_ind + graph_num)"), valign=:bottom, font=:bold)
    resize_to_layout!(figure)
end

function add_param_plots!(figure::Figure, param_plot_data::Tuple, boundary_cond::Symbol)
    J, Δ, ϵ, (J_avg, Δ_avg, ϵ_avg) = param_plot_data

    J_grid = figure[1, 1] = GridLayout()
    Δ_grid = figure[2, 1] = GridLayout()
    ϵ_grid = figure[3, 1] = GridLayout()

    ax_g1 = Axis(J_grid, xlabel = L"\textrm{J}~(J)")
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

    Label(J_grid[Top()], "Fluctuation of the J, Δ, ϵ Parameters", valign=:bottom, font=:bold)

    resize_to_layout!(figure)
end

function bloch_chain_anim()

end

function display_fig(figure::Union{Figure, Vector{Figure}}) # Displaying Makie's images only works in Visual Studio Code
    
    if figure isa Vector
        for fig ∈ figure 
            fig
            readline()
        end
        return nothing
    else
        fig
        readline()
        return nothing
    end
end