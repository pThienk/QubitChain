include("../src/QubitChain.jl")

using .QubitChain, BenchmarkTools

graphs::Vector{BitMatrix} = [[0 0; 0 1], [1 0; 0 1], [1 0; 0 1], [1 0; 0 1], [1 0; 0 1], [1 0; 0 0]]
boundary_cond::Symbol = :open_ended
N::Int64 = 6

println(@btime QubitChain.I_build_DEQsys_expr(N, boundary_cond, graphs))

readline()

graphs_r::Vector{BitMatrix} = [[1 0; 0 1], [1 0; 0 1], [1 0; 0 1], [1 0; 0 1], [1 0; 0 1], [1 0; 0 1]]
boundary_cond_r::Symbol = :round

println(@btime QubitChain.I_build_DEQsys_expr(N, boundary_cond_r, graphs_r))

#= Works as intended =#