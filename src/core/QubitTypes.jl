# Qubit Type System

import Base.@kwdef as @kwdef

abstract type QComponent end

abstract type QOperator end

abstract type QHamiltonian <: QOperator end

@kwdef mutable struct QHamiltonian_S <: QHamiltonian
    J::Float64
    Δ::Float64
    ϵ::Float64
    interaction_graph::BitMatrix = [0 0; 0 0]
end

@kwdef struct QubitData <: QComponent
    dot_1_raw::Vector{ComplexF64}
    dot_2_raw::Vector{ComplexF64}
    dot_1::Vector{ComplexF64}
    dot_2::Vector{ComplexF64}
    t::Vector{Float64}
end

@kwdef mutable struct QChain <: QComponent
    qubits::Vector{QubitData} = []
    N::Int64
    boundary_cond::Symbol = :open_ended
    H::Vector{QHamiltonian_S}
    full_sol = nothing # ODESolution
end

@kwdef struct QChainInitial <: QComponent
    N::Int64
    initial_states::Vector{Tuple{ComplexF64, ComplexF64}}
end

@kwdef struct QChainData <: QComponent
    qubits::Vector{QubitData} = []
    N::Int64
    t_raw::Vector{Float64} = []
    boundary_cond::Symbol = :open_ended
end

@kwdef struct DensityMatrix <: QOperator
    components::Vector{Matrix{ComplexF64}} = []
    t::Vector{Float64}
end
