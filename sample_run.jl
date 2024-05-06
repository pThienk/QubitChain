include("src/QubitChain.jl")

using .QubitChain, Random, Unitful

J0::Float64 = 1.00
Δ0::Float64 = 0.50
ϵ0::Float64 = 0.0

round_num::Int64 = 10

N::Int64 = 20
boundary_cond::Symbol = :open_ended
interaction_model::BitMatrix = [1 0; 0 1]

set_param_y_ticks(0.01)

set_solution_step_size(0.05)

set_t_ticks(20.0)

tspan = (0.0, 200.0)

set_t_unit(u"1e-14 * s")

H = H_rand_model(J0, Δ0, ϵ0, N, boundary_cond, interaction_model; J_fluc=0.0, Δ_fluc=0.0)

initial_cond = initial_chain_model(N, (1:N, ANTI_CAT))

q_chain::QChain = QChain(N = N, boundary_cond = boundary_cond, H = H)
simulate!(q_chain, initial_cond, tspan; solver_settings=(reltol=1e-6,), comf_bypass=true)

visualize(q_chain; save="Sample_Run.png", qubit_plot_settings=(adj_y_lim=true,))

save_chain(q_chain; save="sample_save.csv")