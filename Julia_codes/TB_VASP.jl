using Plots
using DelimitedFiles
using CSV
using Plots.PlotMeasures
using DataFrames
using LinearAlgebra
using LaTeXStrings


const FILE_TB = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/df_2/NbP_hr.dat"
const FILEW_TB = raw"D:\Bandstructure\VASP\NbP\Supercell_final\Coord_1\nSpin_nH\Tight-Binding\df\NbP_tb.dat" # Windows Path

## ## ## Instructions ## ## ##
# Copy and paste the tb.dat file to the directory where DataFrame files will be saved.
# Then, copy the path to the tb.dat file and paste it in the FILE_TB constant.


## Function to read a specific line from a file
function read_line(file_path::String, line_number::Int64)
    open(file_path) do file
        for i in 1:line_number-1
            readline(file)  # Skip the line_number-1 lines
        end
        return readline(file)  # Read the fifth line
    end
end

## Function to load the Tight-Binding Hamiltonian, by reading all supercells at once
function load_file(FILE_PATH_TB::String)
    num_wan = parse(Int64,read_line(FILE_PATH_TB,2))
    num_supcel = parse(Int64,read_line(FILE_PATH_TB,3))
    start_line = Int(4 + ceil(num_supcel/15) )
    df = open(FILE_PATH_TB) do io
        CSV.read(io, DataFrame; 
        header=["Rx","Ry","Rz","n₀", "nᵢ", "En", "iEn"],
        delim=' ', ignorerepeated=true,
        skipto=start_line
        )
    end
    df.iEn[ismissing.(df.iEn)] .= 0
    df.iEn = identity.(df.iEn)
    for i in 1:num_supcel
        CSV.write(FILE_PATH_TB[1:end-4]*"$(i).csv", df[(i-1)*(num_wan^2)+1:(i-1)*(num_wan^2)+num_wan^2,:])
    end
end

# ## Function to load the Tight-Binding Hamiltonian, by reading all supercells at once
# function load_file2(FILE_PATH_TB::String)
#     num_wan = parse(Int64,read_line(FILE_PATH_TB,2))
#     num_supcel = parse(Int64,read_line(FILE_PATH_TB,3))
#     start_line = Int(7 + ceil(num_supcel/15) )
#     df = open(FILE_PATH_TB) do io
#         CSV.read(io, DataFrame; 
#         header=["n₀", "nᵢ", "En", "iEn"],
#         delim=' ', ignorerepeated=true,
#         skipto=start_line,
#         limit=(num_wan^2+1)*num_supcel
#         )
#     end
#     df.iEn[ismissing.(df.iEn)] .= 0
#     df.iEn = identity.(df.iEn)
#     for i in 1:num_supcel
#         CSV.write(FILE_PATH_TB[1:end-4]*"$(i).csv", df[(i-1)*(num_wan^2+1)+1:(i-1)*(num_wan^2+1)+num_wan^2+1,:])
#     end
# end