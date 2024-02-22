using PeriodicTable
using LinearAlgebra
using Statistics
using Plots
using DelimitedFiles
using CSV

# Load data
tst = raw"/home/r_floren/BandStructure/QE/NbP/Local/supercell/tst.in"

function find_species(tst::String)
    s = open(tst) do file
        lines = readlines(file)
        match = findfirst(i -> occursin("ATOMS",i),lines)
        match1 = findfirst(i -> occursin("FRAMES",i),lines)
        full_elmt = map(x->split.(lines[x]),match+1:match1-1)
        return unique(map(x -> parse(Int64,x[1]),full_elmt))
    end
end

function plane_eraser(tst::String,n::Int64)
    file_content = split.(read(tst, String), '\n')
    for i in 1:n
        # Find the number of atoms
        m_line1 = findfirst(i -> occursin("ATOMIC_POSITIONS", i), file_content);
        m_line2 = findfirst(i -> occursin("K_POINTS", i), file_content);
    
        erse = parse.(Float64,map(x->split.(file_content[m_line1+1:m_line2-1])[x][4],1:m_line2-m_line1-1))
        erse = findall(erse .== maximum(erse))
    
        result = setdiff(1:length(file_content),m_line1.+erse)
        file_content = file_content[result];
    end
    
    return file_content
end

tmp = plane_eraser(tst,2)
open(tst, "w") do file
    println(file, join(tmp, '\n'))
end