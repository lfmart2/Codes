using PeriodicTable
using LinearAlgebra
using Statistics
using Plots
using DelimitedFiles
using CSV

# Load data
tst = raw"/home/r_floren/BandStructure/QE/NbP/Local/supercell/tst.in"

function find_min(m::Vector{Float64},n::Int64)
    tmp = copy(m);
    result = [0.0];
    tmp1 = 0.0;
    for i in 1:n
        erse = findall(tmp .== minimum(tmp))
        result1 = zeros(Int64, length(result) + length(erse))
        result1[1:length(result)] .= result
        result1[length(result)+1:end] .= erse
        result = result1
        tmp[result[2:end]] .= trunc(Int,maximum(m))
    end
    return result[2:end]
end

function fix_plane(tst::String,n::Int64)
    s = open(tst) do file
        lines = readlines(file)
        match = findfirst(i -> occursin("ATOMIC_POSITIONS",i),lines)
        match1 = findfirst(i -> occursin("K_POINTS",i),lines)
        tmp = map(last, [split(lines[x]) for x in match + 1:match1 - 1]);
        vec = parse.(Float64, tmp);
        tmp1 = setdiff(match+1:match1-1,find_min(vec,n).+match)
        lines[tmp1] .= lines[tmp1] .* "  1   1   1"
        lines[find_min(vec,n).+match] .= lines[find_min(vec,n).+match] .* "  0   0   0"
        return lines
    end
end
tmp = fix_plane(tst,2)
open(tst, "w") do file
    println(file, join(tmp, '\n'))
end
# unique(map(x -> parse(Int64,x[1]),full_elmt))