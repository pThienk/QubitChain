# File containing the analytical calculations for the dynamic solutions

function cal_density_matrix(qchain::Union{QChain, QChainData}, purity_states::Tuple{Tuple, Tuple})::Vector{DensityMatrix}
    density_matrices::Vector{DensityMatrix} = []
    
    for i ∈ 1:qchain.N

        density_matrix::DensityMatrix = DensityMatrix(t = qchain.qubits[i].t)
        for j ∈ eachindex(qchain.qubits[i].t)
            
            comp11::ComplexF64 = (purity_states[1]⊙(qchain.qubits[i].dot_1[j], qchain.qubits[i].dot_2[j])) * ((qchain.qubits[i].dot_1[j], qchain.qubits[i].dot_2[j])⊙purity_states[1])
            comp12::ComplexF64 = (purity_states[1]⊙(qchain.qubits[i].dot_1[j], qchain.qubits[i].dot_2[j])) * ((qchain.qubits[i].dot_1[j], qchain.qubits[i].dot_2[j])⊙purity_states[2])
            comp21::ComplexF64 = (purity_states[2]⊙(qchain.qubits[i].dot_1[j], qchain.qubits[i].dot_2[j])) * ((qchain.qubits[i].dot_1[j], qchain.qubits[i].dot_2[j])⊙purity_states[1])
            comp22::ComplexF64 = (purity_states[2]⊙(qchain.qubits[i].dot_1[j], qchain.qubits[i].dot_2[j])) * ((qchain.qubits[i].dot_1[j], qchain.qubits[i].dot_2[j])⊙purity_states[2])

            push!(density_matrix.components, [comp11 comp12; comp21 comp22])
            
        end
        
        push!(density_matrices, density_matrix)
    end

    return density_matrices
end

"""

    Note that purity ≡ tr(ρ^2)

"""
function cal_purity(density_matrices::DensityMatrix...)::Vector{Vector{Float64}}

    if isempty(density_matrices)
        throw("The density matrix cannot be null!")
    end

    purities::Vector{Vector{Float64}} = []

    for i ∈ eachindex(density_matrices)

        purity::Vector{Float64} = []
        for j ∈ eachindex(density_matrices[i].t) 
            push!(purity, Float64( tr(density_matrices[i].components[j] ^ 2) ))
    
        end

        push!(purities, purity)
    end
    
    return purities
end

function ⊙(a::Tuple{<:Number, <:Number}, b::Tuple{<:Number, <:Number})
    return conj(a[1])*b[1] + conj(a[2])*b[2]
end

# For other purposes
function ⊙()
    
end