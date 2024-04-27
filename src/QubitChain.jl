module QubitChain

using DifferentialEquations, LaTeXStrings, Unitful, UnitfulLatexify, Latexify, CairoMakie, Printf, SpecialFunctions, Statistics, Distributions, BenchmarkTools,
    RecipesBase, Random, CSV, DataFrames

export H_rand_model, simulate!, save_chain, visualize, initial_chain_model, set_solution_step_size, set_rand_seed, load_data, load_parameters, set_t_ticks,
set_param_y_ticks

export QHamiltonian_S, QubitData, QChain, QChainInitial, QChainData, QComponent, QOperator, QHamiltonian

export CAT, DOT_0, DOT_1, ANTI_CAT, Y_PLUS, Y_MINUS

# Includes types used to hold simulation inputs and store results
include("core/QubitTypes.jl")

# Includes macros for the core system
include("core/CoreMacros.jl")

# Includes the core components of the program: DEQsolver, Qubit Chain Parser, etc.
include("core/Core.jl")

# Includes the CSV file reader/writer components
include("IO/CSVFile.jl")

# Includes interface functions
include("core/Interface.jl")

# Includes the plotting sub-module with plotting functions call
include("visualization/PlotChain.jl")


end # module QubitChain

#= 

    Notes:

    ∘ time-order: 1e-3 (s)
    ∘ initial (0, 1.0) (z-eigenstates) --> results --> initial cond: (correlated), (entangled)
    - Write the save function (Save type: .csv)
    ∘ Write the ploting pipeline (basic) DONE
    
=#
