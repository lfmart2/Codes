using Plots
using DelimitedFiles
using CSV
using PeriodicTable
using LinearAlgebra
using Printf

tst = raw"/home/r_floren/BandStructure/VASP/NbP/POSCAR.tst"

function read_poscar(file::String)
    f = open(file)
    lines = readlines(f)
    close(f)
    lattice = readdlm(IOBuffer(join(lines[3:5], "\n")))
    atoms = split(lines[6])
    atom_types = Int.(readdlm(IOBuffer(lines[7])))
    atom_positions = readdlm(IOBuffer(join(lines[10:end],"\n")))
    return lattice, atoms, atom_types, atom_positions
end

h_atoms = findall(x->x==maximum(read_poscar(tst)[end][:,3]),read_poscar(tst)[end][:,3])
# Read the content of the file and convert to String
file_content = read(tst,String)
lines1 = split(file_content, '\n');
lines1[6] = lines1[6]*" H"
lines1[7] = lines1[7]*"    $(length(h_atoms))"
# Find the line index where "CELL_PARAMETERS" occurs
h_atoms_pos = Array{String}(undef, length(h_atoms))
for i in 1:length(h_atoms)
    h_atoms_pos[i] =@sprintf "    %.16f    %.16f    %.16f    T  T  T     H%03i" read_poscar(tst)[end][h_atoms[i],1] read_poscar(tst)[end][h_atoms[i],2] minimum(read_poscar(tst)[end][:,3])-maximum(diff(sort(read_poscar(tst)[end][:,3])))*0.8 i
end
open(raw"/home/r_floren/BandStructure/VASP/NbP/POSCAR", "w") do file
    println(file, join(lines1,'\n')*join(h_atoms_pos, '\n'))
end 

