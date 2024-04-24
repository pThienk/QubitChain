include("../src/QubitChain.jl")

using .QubitChain, Random, BenchmarkTools

qubits::Vector{QubitData} = [] 
N::Int64 = 20
boundary_cond::Symbol = :open_ended
H::Vector{QHamiltonian_S} = []
interaction_model::BitMatrix = [0 0; 0 0]

t = 0:0.001:1 |> collect
H = H_rand_model(0.25, 0.5, 0.3, N, boundary_cond, interaction_model)

wave_model_1(x) = sin(10π*x)^2
wave_model_2(x) = cos(10π*x)^2

for i ∈ 1:N
    push!(qubits, QubitData(dot_1 = (wave_model_1.(t .+ 0.2 .* rand())), dot_2 = (wave_model_2.(t .- 0.2 .* rand())), dot_1_raw = rand(ComplexF64, 70),
     dot_2_raw = rand(ComplexF64, 70), t = t))
end

q_chain::QChain = QChain(qubits = qubits, N = N, boundary_cond = boundary_cond, H = H)

save_chain(q_chain; save="C:/Users/Porpun/OneDrive/Desktop/QubitChain/save_data/tests/TEST.csv", inc_raw_t=false)