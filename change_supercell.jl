using Plots
using DelimitedFiles
using CSV
using PeriodicTable
using LinearAlgebra
# Add more packages as needed
tst = raw"/home/r_floren/BandStructure/QE/NbP/Local/supercell/tst.xsf"


function rep_uncell(tst::String)
    s = open(tst) do file
        lines = readlines(file)
        match = findfirst(i -> occursin("nunit",i),lines)
        return match isa Nothing ? nothing : map(x->parse(Float64,x),split(lines[match])[2:4])
    end
end
#
function supercell_vec(tst::String)
    s = open(tst) do file
        lines = readlines(file)
        match = findfirst(i -> occursin("CELL_PARAMETERS",i),lines)
        return match isa Nothing ? nothing : hcat(map(yy->map(x->parse(Float64,x),split(lines[match+yy])[1:3]),1:3)...)
    end
end


function matrix_to_string(mat::Matrix{Float64})
    return ["      $(mat[1])      $(mat[1,2])       $(mat[1,3])";"      $(mat[2])      $(mat[2,2])       $(mat[2,3])";"      $(mat[3,1])      $(mat[3,2])       $(mat[3,3])"]
end
# File path
tst1 = raw"/home/r_floren/BandStructure/QE/NbP/Local/supercell/tst.in"

# Read the content of the file and convert to String
file_content = read(tst1, String)

# Find the line index where "CELL_PARAMETERS" occurs
lines = split(file_content, '\n')
m_line = findfirst(i -> occursin("CELL_PARAMETERS", i), lines)

# Check if "CELL_PARAMETERS" is found in the file
if m_line !== nothing
    # Replace the lines from match to match+2 with the new supercell vectors
    lines[m_line] = "CELL_PARAMETERS angstrom"
    lines[m_line+1:m_line+3] .= matrix_to_string(transpose(supercell_vec(tst1)).*(rep_uncell(tst).*[1,1,2]))

    # Open the file in write mode and write the modified content, adding a 
    # \n at the end of every line
    open(tst1, "w") do file
        println(file, join(lines, '\n'))
    end
else
    println("Error: 'CELL_PARAMETERS' not found in the file.")
end
