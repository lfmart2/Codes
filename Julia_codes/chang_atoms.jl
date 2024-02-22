using Plots
using DelimitedFiles
using CSV
using PeriodicTable

tst = raw"/home/r_floren/BandStructure/QE/NbP/Local/supercell/tst.xsf"

function find_species(tst::String)
    s = open(tst) do file
        lines = readlines(file)
        match = findfirst(i -> occursin("ATOMS",i),lines)
        match1 = findfirst(i -> occursin("FRAMES",i),lines)
        full_elmt = map(x->split.(lines[x]),match+1:match1-1)
        return unique(map(x -> parse(Int64,x[1]),full_elmt))
    end
end

function rplc_atoms(tst::String)
    s = open(tst) do file
        lines = readlines(file)
        match = findfirst(i -> occursin("ATOMS",i),lines)
        match1 = findfirst(i -> occursin("FRAMES",i),lines)
        full_elmt = map(x->lines[x],match+1:match1-1)
        for i in 1:length(find_species(tst))
            full_elmt = map(x->replace(full_elmt[x]," $(find_species(tst)[i]) "=>elements[find_species(tst)[i]].symbol),1:length(full_elmt))
        end
        return full_elmt
    end
end

tst1 = raw"/home/r_floren/BandStructure/QE/NbP/Local/supercell/tst.in"

# Read the content of the file and convert to String
file_content = read(tst1, String)

# Find the line index where "CELL_PARAMETERS" occurs
lines1 = split(file_content, '\n');
m_line1 = findfirst(i -> occursin("ATOMIC_POSITIONS", i), lines1);
m_line2 = findfirst(i -> occursin("K_POINTS", i), lines1);
N = length(rplc_atoms(tst)) + length(lines1)-m_line2+m_line1
NN = length(lines1)-1
lines = Array{SubString{String}}(undef, N);
lines[1:m_line1-1] .= lines1[1:m_line1-1];
if m_line1 !== nothing
    # Replace the lines from match to match+2 with the new supercell vectors
    lines[m_line1] = "ATOMIC_POSITIONS angstrom"
    lines[m_line1+1:m_line1+length(rplc_atoms(tst))] .= rplc_atoms(tst)
    lines[m_line1+length(rplc_atoms(tst))+1:N] .= lines1[m_line2:end-1]
    # Open the file in write mode and write the modified content, adding a 
    # \n at the end of every line
    #Replace the correct number of atoms:
    lines[findfirst(i -> occursin("nat = $(m_line2-m_line1-1)", i), lines1)] = "  nat = $(length(rplc_atoms(tst)))"
    open(tst1, "w") do file
        println(file, join(lines, '\n'))
    end

else
    println("Error: 'CELL_PARAMETERS' not found in the file.")
end