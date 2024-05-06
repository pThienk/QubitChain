include("../src/QubitChain.jl")

using .QubitChain, Random, Unitful

J0::Float64 = 1.00
Δ0::Float64 = 0.50
ϵ0::Float64 = 0.0

round_num::Int64 = 10

N::Int64 = 20
boundary_cond::Symbol = :open_ended
interaction_model::BitMatrix = [1 0; 0 1]

set_param_y_ticks(0.01)

set_solution_step_size(0.5)

set_t_ticks(2000.0)

tspan = (0.0, 20000.0)

set_t_unit(u"1e-14 * s")

H = H_rand_model(J0, Δ0, ϵ0, N, boundary_cond, interaction_model; J_fluc=0.0, Δ_fluc=0.0)

initial_cond = initial_chain_model(N, (1:N, ANTI_CAT))

#for i ∈ 1:round_num 
    H[1].J = 0.85 # (rand(MersenneTwister(123456789)) * 2.0 / 5.0) + 0.8

    q_chain::QChain = QChain(N = N, boundary_cond = boundary_cond, H = H)
    simulate!(q_chain, initial_cond, tspan; solver_settings=(reltol=1e-6,), comf_bypass=true)

    visualize(q_chain; save="C:/Users/Porpun/OneDrive/Desktop/QubitChain/visuals/1_FLUC-PARAM_ANTI-CAT_HBAR_REG/1_fluc-param_cat.png", qubit_plot_settings=(adj_y_lim=true,))
#end

plot_purities(q_chain; led_labels=("J", [0.85]), filename="C:/Users/Porpun/OneDrive/Desktop/QubitChain/visuals/1_FLUC-PARAM_ANTI-CAT_HBAR_REG/1_fluc-param_cat.png",
 purities_plot_settings=(adj_y_lims=(0.5, 1.2),))



#save_chain(q_chain; save="sample_save.csv")

# C:/Users/Porpun/OneDrive/Desktop/QubitChain/visuals/1_FLUC-PARAM_CAT_10-SHOT_R-0008-0012/1_fluc-param_cat_$(i).png