include("../src/QubitChain.jl")

using .QubitChain, Random, BenchmarkTools

qubits::Vector{QubitData} = [] 
N::Int64 = 20
boundary_cond::Symbol = :open_ended
H::Vector{QHamiltonian_S} = []

for i ∈ 1:N
    push!(H, QHamiltonian_S(J = (3*rand()), Δ = (3*rand()), ϵ = (3*rand())))
    push!(qubits, QubitData(dot_1 = rand(ComplexF64, 100), dot_2 = rand(ComplexF64, 100), dot_1_raw = rand(ComplexF64, 70),
     dot_2_raw = rand(ComplexF64, 70), t = rand(Float64, 100)))
end

q_chain::QChain = QChain(qubits = qubits, N = N, boundary_cond = boundary_cond, H = H)

@time visualize(q_chain; type=:graph, save="C:/Users/Porpun/OneDrive/Desktop/QubitChain/visuals/test")