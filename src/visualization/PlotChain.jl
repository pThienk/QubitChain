# Collection of plotting interface fucntions and its component calls

# Includes the plot recipes
include("PlotRecipes.jl")

# Includes the animation functions
include("ChainAnim.jl")

function plot_chain(chain_data::Union{QChain, QChainData}; filename::String="", param_plot_settings=(), qubit_plot_settings=())
    qubits_plot_data::Vector{Tuple}, N::Int64, boundary_cond::Symbol = plot_obj_convert(chain_data, :graph)
    acc_plot_num::Int64 = 0

    qubit_figs::Vector{Figure} = []
    param_fig::Figure = Figure()

    with_theme(theme_latexfonts()) do 
        
        param_fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (800, 700))
        
        add_param_plots!(param_fig, param_to_scatter_arr(chain_data), boundary_cond; param_plot_settings...)

        for i ∈ 1:floor(Int64, N/10) 
            push!(qubit_figs, Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size=(800, 1000)))
            add_qubit_plot!(qubit_figs[i], qubits_plot_data, acc_plot_num, boundary_cond; qubit_plot_settings...)
            acc_plot_num += 10
        end

        if N % 10 > 0
            push!(qubit_figs, Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size=(800, 1000)))
            add_qubit_plot!(qubit_figs[i], qubits_plot_data, acc_plot_num, boundary_cond, N % 10; qubit_plot_settings...)
        end

        @info "Finished creating figures. \nDisplaying..."

        # display_fig(param_fig)
        # display_fig(qubit_figs)

    end

    if filename ≠ ""
        save(filename * "_param.png", param_fig)

        for i ∈ eachindex(qubit_figs)
            save(filename * "_figure_" * string(('a' + i - 1)) * ".png", qubit_figs[i])
        end

        @info "Saving completed."
    end
    
end

function add_qubit_plot!(figure::Figure, plot_data::Vector{Tuple}, plot_ind::Int64, boundary_cond::Symbol, graph_num::Int64=10; settings...)
    
    w_grid = figure[1, 1] = GridLayout()
    for i ∈ 1:graph_num

        t::Vector{Float64}, dot_1_data::Vector{Float64}, dot_2_data::Vector{Float64} = plot_data[i + plot_ind]

        if :disp_inv ∈ keys(settings)
            δt = t[2] - t[1]

            t = t[(floor(Int64, settings[:disp_inv][1]/δt) + 1):(floor(Int64, settings[:disp_inv][2]/δt) + 1)]
            dot_1_data = dot_1_data[(floor(Int64, settings[:disp_inv][1]/δt) + 1):(floor(Int64, settings[:disp_inv][2]/δt) + 1)]
            dot_2_data = dot_2_data[(floor(Int64, settings[:disp_inv][1]/δt) + 1):(floor(Int64, settings[:disp_inv][2]/δt) + 1)]
        end

        grid_lvl = w_grid[i, 1:2] = GridLayout()
        grid_left = grid_lvl[1, 1] = GridLayout()
        grid_right = grid_lvl[1, 2] = GridLayout()

        ax_left = Axis(grid_left[1, 1], xlabel = L"t ($\times$ %$(latexify(time_unit)))", ylabel = "qubit $(plot_ind + i)")

        if :plot_type ∈ keys(settings) && settings[:plot_type] == :only_dot0
            lines!(ax_left, t, dot_1_data, color=:blue, label="dot 0")
        elseif :plot_type ∈ keys(settings) && settings[:plot_type] == :only_dot1
            lines!(ax_left, t, dot_2_data, color=:red, label="dot 1")
        else
            lines!(ax_left, t, dot_1_data, color=:blue, label="dot 0")
            lines!(ax_left, t, dot_2_data, color=:red, label="dot 1")
        end

        if i ≠ graph_num
            hidexdecorations!(ax_left, grid = false)
        end

        if :adj_y_lim ∈ keys(settings) && settings[:adj_y_lim] == true
            xlims!(ax_left, low = t[1])
            ylims!(ax_left, low = min([dot_1_data; dot_2_data]...))
        else
            xlims!(ax_left, low = t[1])
            ylims!(ax_left, low = 0)
        end

        ax_left.xticks = t[1]:t_ticks:(max(t...) + t_ticks)

        Legend(grid_right[1, 1], ax_left)
        
    end

    colgap!(w_grid, 10)
    rowgap!(w_grid, 15)

    Label(w_grid[1, 1:2, Top()], "Dynamic of qubit $(plot_ind + 1)-$(plot_ind + graph_num) with boundary condition: $(string(boundary_cond))",
     valign=:bottom, halign=:center, font=:bold, padding=(0, 0, 5, 0))

    resize_to_layout!(figure)
end

function add_param_plots!(figure::Figure, param_plot_data::Tuple, boundary_cond::Symbol; settings...)
    J, Δ, ϵ, (J_avg, Δ_avg, ϵ_avg) = param_plot_data

    grid = figure[1, 1] = GridLayout()

    J_grid = grid[1, 1] = GridLayout()
    Δ_grid = grid[2, 1] = GridLayout()
    ϵ_grid = grid[3, 1] = GridLayout()

    ax_g1 = Axis(J_grid[1, 1], ylabel = "J (Hz in ħ unit)", ytickformat = "{:.4f}")
    ax_g2 = Axis(Δ_grid[1, 1], ylabel = "Δ (Hz in ħ unit)", ytickformat = "{:.4f}")
    ax_g3 = Axis(ϵ_grid[1, 1], ylabel = "ϵ (Hz in ħ unit)", xlabel = L"N", ytickformat = "{:.4f}")

    if boundary_cond == :open_ended
        J_points = Point2f.((1:length(J) |> collect)[1:end-1], J[1:end-1]) 
        scatter!(ax_g1, J_points , label="J", color=:purple)
    else
        J_points = Point2f.((1:length(J) |> collect), J) 
        scatter!(ax_g1, J_points, label="J", color=:purple)
    end

    Δ_points = Point2f.((1:length(Δ) |> collect), Δ)
    ϵ_points = Point2f.((1:length(ϵ) |> collect), ϵ)
    
    scatter!(ax_g2, Δ_points, label="Δ", color=:blue)
    scatter!(ax_g3, ϵ_points, label="ϵ", color=:gold)

    lines!(ax_g1, 0..length(J), x -> J_avg, linestyle=:dash, color=:black, alpha=0.75)
    lines!(ax_g2, 0..length(Δ), x -> Δ_avg, linestyle=:dash, color=:black, alpha=0.75)
    lines!(ax_g3, 0..length(ϵ), x -> ϵ_avg, linestyle=:dash, color=:black, alpha=0.75)

    ax_g1.yticks = min(J...):param_y_ticks:max(J...)
    ax_g2.yticks = min(Δ...):param_y_ticks:max(Δ...)
    ax_g3.yticks = min(ϵ...):param_y_ticks:max(ϵ...)

    ax_g1.xticks = 1:length(J)
    ax_g2.xticks = 1:length(Δ)
    ax_g3.xticks = 1:length(ϵ)

    hidexdecorations!(ax_g1, grid = false)
    hidexdecorations!(ax_g2, grid = false)

    xlims!(ax_g1, low = 0, high = length(J) + 1)
    xlims!(ax_g2, low = 0, high = length(Δ) + 1)
    xlims!(ax_g3, low = 0, high = length(ϵ) + 1)

    rowgap!(grid, 7)

    Label(grid[1, 1, Top()], "Fluctuation of the J, Δ, ϵ Parameters with boundary condition: $(string(boundary_cond))",
     valign=:bottom, halign=:center, font=:bold, padding=(0, 0, 5, 0))

    resize_to_layout!(figure)
end

function display_fig(figure::Union{Figure, Vector{Figure}}) # Displaying Makie's images only works in Visual Studio Code
    
    if figure isa Vector
        for fig ∈ figure 
            fig
        end
    else
        figure
    end
end