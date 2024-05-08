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

function cal_sub_density_matrix_at_t(qchain::Union{QChain, QChainData}, purity_states::Tuple{Tuple, Tuple}, indexes::Tuple{Int64, Int64, Int64})::Matrix{ComplexF64}

    i::Int64, j::Int64, k::Int64 = indexes

    sub_density_matrix_t::Matrix{ComplexF64} = zeros(2, 2)

    comp11::ComplexF64 = (purity_states[1]⊙(qchain.qubits[i].dot_1[k], qchain.qubits[i].dot_2[k])) * ((qchain.qubits[j].dot_1[k], qchain.qubits[i].dot_2[k])⊙purity_states[1])
    comp12::ComplexF64 = (purity_states[1]⊙(qchain.qubits[i].dot_1[k], qchain.qubits[i].dot_2[k])) * ((qchain.qubits[j].dot_1[k], qchain.qubits[i].dot_2[k])⊙purity_states[2])
    comp21::ComplexF64 = (purity_states[2]⊙(qchain.qubits[i].dot_1[k], qchain.qubits[i].dot_2[k])) * ((qchain.qubits[j].dot_1[k], qchain.qubits[i].dot_2[k])⊙purity_states[1])
    comp22::ComplexF64 = (purity_states[2]⊙(qchain.qubits[i].dot_1[k], qchain.qubits[i].dot_2[k])) * ((qchain.qubits[j].dot_1[k], qchain.qubits[i].dot_2[k])⊙purity_states[2])

    sub_density_matrix_t =  [comp11 comp12; comp21 comp22]

    return sub_density_matrix_t
end

function cal_system_purity(qchain::Union{QChain, QChainData}, purity_states::Tuple{Tuple, Tuple})::Tuple{Vector{Float64}, Vector{Float64}}

    system_density_matrix::DensityMatrix = DensityMatrix(t = qchain.qubits[1].t)

    for k ∈ eachindex(qchain.qubits[1].t)

        system_density_matrix_rows::Vector{Matrix{ComplexF64}} = []
        for i ∈ 1:qchain.N

            row_density_matrix_comp::Vector{Matrix{ComplexF64}} = []
            for j ∈ 1:qchain.N 

                push!(row_density_matrix_comp, cal_sub_density_matrix_at_t(qchain, purity_states, (i, j, k)))
            end

            push!(system_density_matrix_rows, hcat(row_density_matrix_comp...))
        end

        push!(system_density_matrix.components, vcat(system_density_matrix_rows...))
    end

    purities::Vector{Float64} = real.( tr.( system_density_matrix.components .^ 2 ) )

    return (purities, system_density_matrix.t)
end

"""

    Note that purity ≡ tr(ρ^2)

"""
function cal_purity(density_matrices::Vector{DensityMatrix}...)::Vector{Vector{Tuple{Vector{Float64}, Vector{Float64}}}}

    if isempty(density_matrices)
        throw("The density matrix cannot be null!")
    end

    purities::Vector{Vector{Tuple{Vector{Float64}, Vector{Float64}}}} = []

    for i ∈ eachindex(density_matrices[1])

        chain_purities::Vector{Tuple{Vector{Float64}, Vector{Float64}}} = []
        for j ∈ eachindex(density_matrices)
           
            purity::Vector{Float64} = real.( tr.(density_matrices[j][i].components .^ 2) )
            
            push!(chain_purities, (purity, density_matrices[j][i].t))

        end

        push!(purities, chain_purities)
    end
    
    return purities
end

function ⊙(a::Tuple{<:Number, <:Number}, b::Tuple{<:Number, <:Number})
    return conj(a[1])*b[1] + conj(a[2])*b[2]
end

# For other purposes
function ⊙()
    
end