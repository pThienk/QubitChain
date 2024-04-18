# Contains function for handling writing and reading .csv file

# Includes abstract file writing/reading
include("FileIO.jl")

function write_to_csv(filename::String, header::Vector{String}, data::Vector...)

    if isempty(data)
        throw("Invalid: data must not be empty!")
    end

    
    row_array::Vector{String} = []

    push!(row_array, join(header, ","))
    # println(max((length.(data))...))
    for i ∈ 1:max((length.(data))...)
        row_str::String = ""
        for j ∈ eachindex(header) 
            if length(data[j]) < i
                if j == eachindex(header)[end]
                    row_str = join([row_str, ""])
                else
                    row_str = join([row_str, ","])
                end
            else # Fix comma problem
                if j == eachindex(header)[end]
                    row_str = join([row_str, string(data[j][i])])
                else
                    row_str = join([row_str, string(data[j][i]), ","])
                end 
            end
        end
        
        push!(row_array, row_str)
    end

    write_file(filename, row_array...)
end

function read_from_csv(filename::String)::DataFrame
    return CSV.File(filename; normalizenames=true) |> DataFrame
end