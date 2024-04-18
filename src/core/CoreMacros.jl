# Macros and macro wrappers for supplementing code in the core

function build_DEQ_sys(N::Int64, boundary_cond::Symbol, interaction_graphs::Vector{BitMatrix})
    return eval(Meta.parse(I_build_DEQsys_expr(N, boundary_cond, interaction_graphs)))
 end

 macro logger()
    # Logging system
 end

"""
    ***Internal Function: Not for end user***
"""
 function I_build_DEQsys_expr(N::Int64, boundary_cond::Symbol, graphs::Vector{BitMatrix})::String
    param_str::String = ""
    comp_str::String = ""

    deq_array::Vector{String} = []
    deq_ind::Int64 = 1
    for i ∈ 1:N
        if i ≠ N
            param_str = join([param_str, "J", i, "_", i+1, ",",  " ", "Δ", i, ",", " ", "ϵ", i, ",", " "])
        elseif boundary_cond == :open_ended
            param_str = join([param_str, "J", i, "_", i+1, ",",  " ", "Δ", i, ",", " ", "ϵ", i, ",", " "])
        elseif boundary_cond == :round
            param_str = join([param_str, "J", i, "_", 1, ",",  " ", "Δ", i, ",", " ", "ϵ", i, ",", " "])
        end

        deq_temp_str_arr::Vector{String} = (join(["ϵ", i])⊗parse_str_matrix(σz) ⊕ 
            join(["Δ", i])⊗parse_str_matrix(σx) ⊕ mask_interaction_param_S(i, graphs[i], N))⊗[join(["A", i, "_", 0]), join(["A", i, "_", 1])]

        push!(deq_array, join(["du[", deq_ind, "] = ", "(-im/ħ)*(", deq_temp_str_arr[1], ")"]), 
            join(["du[", deq_ind+1, "] = ", "(-im/ħ)*(", deq_temp_str_arr[2], ")"]))
        
        deq_ind += 2
            
        comp_str = join([comp_str, "A", i, "_", 0, ",", " ", "A", i, "_", 1, ",", " "])
    end

    param_str *= "= p"
    comp_str *= "= u"

    deq_sys_str::String = join(deq_array, "\n")

    return join(["function(du, u, p, t)", "\n", param_str, "\n", comp_str, "\n", deq_sys_str, "\n", "end"])
 end
