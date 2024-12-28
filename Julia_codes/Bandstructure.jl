using Plots
using DelimitedFiles
using HDF5
using CSV
using Plots.PlotMeasures
using DataFrames
using LinearAlgebra
using LaTeXStrings


# The following code calculates the bandstructure of a material using the PROCAR file from VASP. In this code, 
# the bandstructure is plotted with the energy levels of the bands and the Fermi energy. The code also includes
# the density of states projected on the d orbitals of the material. The code also includes the symmetry points.
######## Details of the PROCAR file ########
# The PROCAR file contains the band structure of the material. The file contains the k-points data[1,2], 
# the energy levels of the bands data[2,5], and the total and projected density of states of all atoms sum(data[12,6:10]).
# The PROCAR file is divided into two blocks, the first block contains the data of the spin up electrons and the second block
# contains the data of the spin down electrons. The data is divided into blocks of size n_bands+1, where n_bands is the number of bands.

# The following code reads the PROCAR file, the DOSCAR file, and the KPOINTS file for the bandstructure calculation of
# non-polarized withouth Hydrogen atom.
const FILE_PROCAR_nSpin_nH1    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Bandstructure/Bandstructure1/PROCAR"
const FILE_PROCAR_nSpin_nH2    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Bandstructure/Bandstructure2/PROCAR"
const FILE_DOSCAR_nSpin_nH1    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Bandstructure/Bandstructure1/DOSCAR"
const FILE_DOSCAR_nSpin_nH2    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Bandstructure/Bandstructure2/DOSCAR"
const PATH_DF_nSpin_nH         = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Bandstructure/df/"
const FILE_BANDS_INFO_nSpin_nH = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Bandstructure/df/KPOINTS.tot"

# The following code reads the PROCAR file, the DOSCAR file, and the KPOINTS file for the bandstructure calculation of
# non-polarized with Hydrogen atom.
const FILE_PROCAR_nSpin_H1    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_H/Bandstructure/Bandstructure1/PROCAR"
const FILE_PROCAR_nSpin_H2    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_H/Bandstructure/Bandstructure2/PROCAR"
const FILE_DOSCAR_nSpin_H1    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_H/Bandstructure/Bandstructure1/DOSCAR"
const FILE_DOSCAR_nSpin_H2    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_H/Bandstructure/Bandstructure2/DOSCAR"
const PATH_DF_nSpin_H         = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_H/Bandstructure/df/"
const FILE_BANDS_INFO_nSpin_H = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_H/Bandstructure/df/KPOINTS.tot"

# The following code reads the PROCAR file, the DOSCAR file, and the KPOINTS file for the bandstructure calculation of
# spin-polarized withouth Hydrogen atom.
const FILE_PROCAR_Spin_nH1    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/Bandstructure1/PROCAR"
const FILE_PROCAR_Spin_nH2    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/Bandstructure2/PROCAR"
const FILE_PROCAR_Spin_nH3    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/Bandstructure3/PROCAR"
const FILE_PROCAR_Spin_nH4    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/Bandstructure4/PROCAR"
const FILE_DOSCAR_Spin_nH1    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/Bandstructure1/DOSCAR"
const FILE_DOSCAR_Spin_nH2    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/Bandstructure2/DOSCAR"
const FILE_DOSCAR_Spin_nH3    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/Bandstructure3/DOSCAR"
const FILE_DOSCAR_Spin_nH4    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/Bandstructure4/DOSCAR"
const PATH_DF_Spin_nH         = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/df/"
const FILE_BANDS_INFO_Spin_nH = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/df/KPOINTS.tot"

# The following code reads the PROCAR file, the DOSCAR file, and the KPOINTS file for the bandstructure calculation of
# spin-polarized with Hydrogen atom.
const FILE_PROCAR_Spin_H1    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/Bandstructure1/PROCAR"
const FILE_PROCAR_Spin_H2    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/Bandstructure2/PROCAR"
const FILE_PROCAR_Spin_H3    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/Bandstructure3/PROCAR"
const FILE_PROCAR_Spin_H4    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/Bandstructure4/PROCAR"
const FILE_DOSCAR_Spin_H1    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/Bandstructure1/DOSCAR"
const FILE_DOSCAR_Spin_H2    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/Bandstructure2/DOSCAR"
const FILE_DOSCAR_Spin_H3    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/Bandstructure3/DOSCAR"
const FILE_DOSCAR_Spin_H4    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/Bandstructure4/DOSCAR"
const PATH_DF_Spin_H         = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/df/"
const FILE_BANDS_INFO_Spin_H = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/df/KPOINTS.tot"


################################################################################################################
############################################## CODE STARTS HERE ################################################
################################################################################################################

# The following code reads the PROCAR file, the DOSCAR file for the creation of the DataFrames with the bandstructure.
# This process is fundamental, since we need to "chunk" the Bandstructure Calculations (VASP) for the memory usage.
# Thus, we first need to substract the information of the PROCAR file and the Fermi energy from the DOSCAR file, and 
# then we can create the DataFrames. The DataFrames are saved in the df folder. The format that the DataFrames are saved
# already includes the projections of the atomic orbitals and the atoms. The DataFrames are saved in the .csv format.
# An example to run this is the following: filter_file(FILE_PROCAR_nSpin_nH1, FILE_DOSCAR_nSpin_nH1)
function filter_file(FILE_PATH_BAND::String, FILE_PATH_FE::String)
    # Open the PROCAR file and read it into a DataFrame, skipping the first row and ignoring empty rows
    df = open(FILE_PATH_BAND) do io
        dropmissing(CSV.read(io, DataFrame; header=["ion", "s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2", "tot"], skipto =2, ignoreemptyrows=true, delim=' ', ignorerepeated=true, footerskip=1), :ion)
    end

    # Remove rows where the "ion" column has values "k-point" or "ion"
    df = deleteat!(df, findall(i -> (i == "k-point") || (i == "ion"), df.ion))

    # Check if there is only one "#" in the "ion" column
    if length(findall(i -> i == "#", df.ion)) == 1
        # Extract initial information from the DataFrame
        init_info = findall(i -> i == "#", df.ion)[1]
        info = DataFrame(band = ["Non-polarized", "#-kpoints:", "#-bands:", "#-ions:", "E-Fermi:"],
            energy = ["", parse(Int64, df[init_info, 4]), parse(Int64, df[init_info, 8]), df[init_info, 12], extract_fermi_energy(FILE_PATH_FE)])

        # Parse the energies from the DataFrame
        energies = parse.(Float64, df.px[findall(i -> i == "band", df.ion)])
        tmp = DataFrame("band" => repeat(1:info.energy[3], info.energy[2]), "energy" => energies)

        # Remove rows where the "ion" column has values "#" or "band"
        df = select(deleteat!(df, findall(i -> i == "#" || i == "band", df.ion)), collect(1:11))

        # Replace "tot" values in the "ion" column with "0"
        df.ion[findall(i -> i == "tot", df.ion)] .= "0"

        # Parse the columns to appropriate data types
        col = names(df)
        df[!, col[1]] .= parse.(Int64, df[!, col[1]])
        for i in 2:length(col)
            df[!, col[i]] .= parse.(Float64, df[!, col[i]])
        end

        # Create a new directory and folder df:
        n_d = FILE_PATH_BAND[1:findall(isequal('/'), FILE_PATH_BAND)[end-1]] * "df/"
        
        # Write the DataFrame to a .csv file with the name label/tag of the POSCAR file
        CSV.write(n_d * "df" * FILE_PATH_BAND[end-7] * "_bnd.csv", df)
        CSV.write(n_d * "df" * FILE_PATH_BAND[end-7] * "_energies.csv", vcat(info, tmp))
    else

        # Extract information for the positions from where starts/ends each spin polarization
        spn_up_info = findall(i -> i == "#", df.ion)[1]
        spn_dn_info = findall(i -> i == "#", df.ion)[2]

        n_bnds   = 552
        k_points =  parse(Int64, df[spn_up_info, 4])
        n_ions   = df[spn_dn_info, 12]
        
        # Create DataFrames for spin-up and spin-down information
        # info_up = DataFrame(band = ["S-up:", "#-kpoints:", "#-bands:", "#-ions:", "E-Fermi:"],
        #     energy = ["", parse(Int64, df[spn_up_info, 4]), parse(Int64, df[spn_up_info, 8]), df[spn_up_info, 12], extract_fermi_energy(FILE_PATH_FE)])
        # info_dn = DataFrame(band = ["S-dn:", "#-kpoints:", "#-bands:", "#-ions:", "E-Fermi:"],
        #     energy = ["", parse(Int64, df[spn_dn_info, 4]), parse(Int64, df[spn_dn_info, 8]), df[spn_dn_info, 12], extract_fermi_energy(FILE_PATH_FE)])
        info_up = DataFrame(band = ["S-up:", "#-kpoints:", "#-bands:", "#-ions:", "E-Fermi:"],
            energy = ["",k_points, n_bnds, n_ions, extract_fermi_energy(FILE_PATH_FE)])
        info_dn = DataFrame(band = ["S-dn:", "#-kpoints:", "#-bands:", "#-ions:", "E-Fermi:"],
            energy = ["", k_points, n_bnds, n_ions, extract_fermi_energy(FILE_PATH_FE)])
        # Parse the energies for spin-up and spin-down
        # energies_up = parse.(Float64, df.px[findall(i -> i == "band", df.ion)])[1:end÷2]
        # energies_dn = parse.(Float64, df.px[findall(i -> i == "band", df.ion)])[end÷2+1:end]
        energies_up = parse.(Float64, df.px[findall(i -> i == "band", df.ion)])[1:n_bnds*k_points]
        energies_dn = parse.(Float64, df.px[findall(i -> i == "band", df.ion)])[end÷2+1:end÷2+n_bnds*k_points]
        
        # Create temporary DataFrames for spin-up and spin-down energies
        tmp_up = DataFrame("band" => repeat(1:info_up.energy[3], info_up.energy[2]), "energy" => energies_up)
        tmp_dn = DataFrame("band" => repeat(1:info_dn.energy[3], info_dn.energy[2]), "energy" => energies_dn)
        
        # Remove rows where the "ion" column has values "#" or "band"
        df = select(deleteat!(df, findall(i -> i == "#" || i == "band", df.ion)), collect(1:11))
        
        # Replace "tot" values in the "ion" column with "0"
        df.ion[findall(i -> i == "tot", df.ion)] .= "0"
        
        # Parse the columns to appropriate data types
        col = names(df)
        df[!, col[1]] .= parse.(Int64, df[!, col[1]])
        for i in 2:length(col)
            df[!, col[i]] .= parse.(Float64, df[!, col[i]])
        end
        
        # Create a new directory and folder df:
        n_d = FILE_PATH_BAND[1:findall(isequal('/'), FILE_PATH_BAND)[end-1]] * "df/"
        
        # Write the DataFrame to .csv files for spin-up and spin-down
        # CSV.write(n_d * "df" * FILE_PATH_BAND[end-7] * "_bnd_up.csv", df[1:end÷2, :])
        # CSV.write(n_d * "df" * FILE_PATH_BAND[end-7] * "_bnd_down.csv", df[end÷2+1:end, :])
        CSV.write(n_d * "df" * FILE_PATH_BAND[end-7] * "_bnd_up.csv", df[1:(n_ions+1)*n_bnds*k_points, :])
        CSV.write(n_d * "df" * FILE_PATH_BAND[end-7] * "_bnd_down.csv", df[end÷2+1:end÷2+(n_ions+1)*n_bnds*k_points, :])
        
        # Write the energies DataFrames to .csv files for spin-up and spin-down
        CSV.write(n_d * "df" * FILE_PATH_BAND[end-7] * "_energies_up.csv", vcat(info_up, tmp_up))
        CSV.write(n_d * "df" * FILE_PATH_BAND[end-7] * "_energies_down.csv", vcat(info_dn, tmp_dn))
        end
end

# Fucntion to merge two Bandstructure plots, using a PROCAR file with the spin-up and spin-down data separated.
# Function to concatenate the PROCAR files for non-spin-polarized calculations
function cat_DF(DF_PATH::String)
    reader = readdir(DF_PATH)
    if any(endswith.(reader, "up.csv"))
        ## Merging the spin-up Bandstructure data
        tmp = DataFrame()
        for i in reader[endswith.(reader,"bnd_up.csv")]
            tmp1 = open(DF_PATH*i) do io
                CSV.read(io, DataFrame)
            end
            tmp = vcat(tmp, tmp1)
        end
        CSV.write(DF_PATH*"df_bands_total_up.csv", tmp)

        ## Merging the spin-down Bandstructure data
        tmp = DataFrame()
        for i in reader[endswith.(reader,"bnd_down.csv")]
            tmp1 = open(DF_PATH*i) do io
                CSV.read(io, DataFrame)
            end
            tmp = vcat(tmp, tmp1)
        end
        CSV.write(DF_PATH*"df_bands_total_down.csv", tmp)

        # Merging spin-up Energies
        tmp = DataFrame()
        kpoints = 0
        header = DataFrame()
        for i in reader[endswith.(reader,"energies_up.csv")]
            kpoints += open(DF_PATH*i) do io
                CSV.read(io, DataFrame; limit=5).energy[2]
            end
            tmp1 = open(DF_PATH*i) do io
                CSV.read(io, DataFrame; skipto = 7)
            end
            tmp = vcat(tmp, tmp1)
        end
        header = open(DF_PATH*"df1_energies_up.csv") do io
            CSV.read(io, DataFrame; limit=5)
        end
        header.energy[2] = kpoints
        # Writing the energies
        CSV.write(DF_PATH*"df_energies_total_up.csv", vcat(header,tmp))

        # Merging spin-down Energies
        tmp = DataFrame()
        kpoints = 0
        header = DataFrame()
        for i in reader[endswith.(reader,"energies_down.csv")]
            kpoints += open(DF_PATH*i) do io
                CSV.read(io, DataFrame; limit=5).energy[2]
            end
            tmp1 = open(DF_PATH*i) do io
                CSV.read(io, DataFrame; skipto = 7)
            end
            tmp = vcat(tmp, tmp1)
        end
        header = open(DF_PATH*"df1_energies_down.csv") do io
            CSV.read(io, DataFrame; limit=5)
        end
        header.energy[2] = kpoints
        # Writing the energies
        CSV.write(DF_PATH*"df_energies_total_down.csv", vcat(header,tmp))
    else
        tmp = DataFrame()
        header = DataFrame()
        kpoints = 0

        ## Merging the Bandstructure data
        for i in reader[endswith.(reader,"bnd.csv")]
            tmp1 = open(DF_PATH*i) do io
                CSV.read(io, DataFrame)
            end
            tmp = vcat(tmp, tmp1)
        end
        CSV.write(DF_PATH*"df_bands_total.csv", tmp)

        # Merging the Energies
        tmp = DataFrame()
        for i in reader[endswith.(reader,"energies.csv")]
            kpoints += open(DF_PATH*i) do io
                CSV.read(io, DataFrame; limit=5).energy[2]
            end
            tmp1 = open(DF_PATH*i) do io
                CSV.read(io, DataFrame; skipto = 7)
            end
            tmp = vcat(tmp, tmp1)
        end
        header = open(DF_PATH*"df1_energies.csv") do io
            CSV.read(io, DataFrame; limit=5)
        end
        header.energy[2] = kpoints

        # Writing the energies
        CSV.write(DF_PATH*"df_energies_total.csv", vcat(header,tmp))
    end
end

# Code to extract the Fermi energy from the DOSCAR file. The input is the path to the DOSCAR file, and the output is the Fermi energy.
function extract_fermi_energy(FILE_PATH::String)::Float64
    open(FILE_PATH, "r") do file
        for i in 1:5
            readline(file)  # Skip the first 5 lines
        end
        line = readline(file)  # Read the 6th line
        return parse(Float64,split(line)[4])  # Print the 6th line
    end 
end

# Code to extract the symmetry point labels from the KPOINTS file. The input is the path to the KPOINTS file, and the 
# output is an array with the labels of the symmetry points.
function label_sym_points(FILE_PATH::String)
    # Read the symmetry point labels from the KPOINTS file, starting from the 5th row and taking the last column
    data = readdlm(FILE_PATH)[5:end, end]
    
    # Initialize an empty array to store the formatted labels
    x_label = String[]
    
    # Check if the first label is "Gamma" and format it accordingly
    if data[1] == "Gamma"
        x_label = push!(x_label, "\$ \\" * data[1] * "\$")
    else
        x_label = push!(x_label, "\$" * data[1] * "\$")
    end
    # Loop through the symmetry points data in pairs
    for i=2:2:length(data)-1
        if data[i] == data[i+1]
            # If the current and next points are the same
            if data[i] == "Gamma"
                # Format "Gamma" with LaTeX notation
                x_label = push!(x_label,"\$ \\"* data[i]* "\$")
            else
                # Format other points with LaTeX notation
                x_label = push!(x_label,"\$"* data[i]* "\$")
            end
        else
            # If the current and next points are different
            if data[i] == "Gamma"
                # Format "Gamma" with LaTeX notation and add a separator
                x_label = push!(x_label,"\$ \\"* data[i]* "|"* data[i+1]* "\$")
            elseif data[i+1] == "Gamma"
                # Format the next point as "Gamma" with LaTeX notation and add a separator
                x_label = push!(x_label,"\$"* data[i]* "|\\"* data[i+1]* "\$")
            else
                # Format both points with LaTeX notation and add a separator
                x_label = push!(x_label,"\$"* data[i]* "|"* data[i+1]* "\$")
            end
        end
    end

    # Check and format the last point
    if data[end] == "Gamma"
        x_label = push!(x_label,"\$ \\"* data[end]* "\$")
    else
        x_label = push!(x_label,"\$"* data[end]* "\$")
    end

end

# Code to explicity select the atomic orbitals to be plotted. The input is a vector with the desired atomic orbitals.
# For instance, the input "p" will give and output of ["ion","px","py","pz"]. The output is an array with the atomic orbitals.
function atomic_orbital(a_orb::Vector{String})
    a_orbital = String[]
    a_orbital = push!(a_orbital,"ion")
    for i in a_orb
        if i == "s"
            a_orbital = push!(a_orbital,"s")
        elseif i == "p"
            a_orbital = push!(a_orbital,"px","py","pz")
        elseif i == "d"
            a_orbital = push!(a_orbital,"dxy","dyz","dz2","dxz","dx2")
        end
    end
    return a_orbital
end


# Function to load and select the specific atomic orbitals and atoms from the PROCAR.csv file. The input is the path to the PROCAR.csv file,
function load_filtered_data(FILE_PATH_BS, FILE_PATH_ENER, select_orbital, select_atom)
    n_kpoints = Int(CSV.File(open(FILE_PATH_ENER); limit=5).energy[2])
    n_bnds = Int(CSV.File(open(FILE_PATH_ENER); limit=5).energy[3])
    # n_atoms = (CSV.File(open(FILE_PATH_ENER); limit=1, header=false, delim=" ", ignorerepeated = true, skipto=2).Column12)[1]
    reader = open(FILE_PATH_BS) do io
        CSV.read(io, DataFrame; select = atomic_orbital(select_orbital))
    end
    return n_kpoints, n_bnds, CSV.File(open(FILE_PATH_ENER); skipto=7).energy, select!(filter(AsTable(:) => nt -> nt.ion in select_atom, reader), Not([:ion]))
end

# Code to plot the bandstructure of the material for spin-up and spin-down. The Bandstructure lines are colored (opacity)
# according to the density of states projected on the chosen orbital of the chosen atom.
# The input files are the paths to the PROCAR file, the KPOINTS file (FILE_PATH_SYM). The next inputs are a vector with the desired atoms to be plotted,
# a vector with the desired orbitals to be plotted, and a vector with the desired spins to be plotted.
# For instance, if I want to plot the bandstructure of the d orbitals of the first and second atoms for the spin-up and 
# s-p orbitals of the second atom for the spin-down, the inputs are n_atom = [[1,2],[2]], a_orbital = [["d"],["s","p"]], and spin = ["up","down"].
# Note that the number labels of the n_atom are according with the POSCAR file, i.e., the first half of atoms are Nb, and the second half are P.
# The output are two plots with the weighted Projected Bandstructure of the material. The plot includes the Fermi energy and the symmetry points.
################### Parameters to tune the plot ###################
# LINE_WIDTH: Array with the width of the lines of the plot. The first element is the width of the bandstructure lines, 
# the second element is the width of the Fermi energy line, and the third element is the vertical symmetry lines.
# LINE_COLOR: Array with the color of the lines of the plot. The first element is the color of the spin-up bandstructure lines,
# the second element is the color of the spin-down bandstructure lines, and the third element is the color of Fermi line.

function Bandstructure_spin_plot(PATH_DF::String,FILE_PATH_SYM::String,n_atom::Vector{Vector{Int64}},spin::Vector{String},a_orbital::Vector{Vector{String}})
    p1 = plot();
        # Parameters for the plot
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:red, :blue, :black]
    ANNOTATE_POSITION = 0.05
    ANNOTATE_SIZE = 9

    fermi_energy = CSV.File(open(PATH_DF*"df_energies_total_up.csv"); limit=5).energy[5]
    sym_points = readdlm(FILE_PATH_SYM)[5:end,1:3]
    
    # First bandstructure
    n_kpoints, n_bnds, energies_1, data_1 = load_filtered_data(
        PATH_DF*"df_bands_total_$(spin[1]).csv",
        PATH_DF*"df_energies_total_$(spin[1]).csv", 
        a_orbital[1], n_atom[1]
        )
    tmp = sum(eachcol(data_1))
    tmp_reshape = reshape(tmp, length(n_atom[1]),:)
    tmp = vec(sum(tmp_reshape, dims=1))
    bands_1 = reshape(tmp, n_bnds, n_kpoints)'
    # Second bandstructure
    n_kpoints, n_bnds, energies_2, data_2 = load_filtered_data(
        PATH_DF*"df_bands_total_$(spin[2]).csv",
        PATH_DF*"df_energies_total_$(spin[2]).csv",
        a_orbital[2], n_atom[2]
        )
    tmp = sum(eachcol(data_2))
    tmp_reshape = reshape(tmp, length(n_atom[2]),:)
    tmp = vec(sum(tmp_reshape, dims=1))
    bands_2 = reshape(tmp, n_bnds, n_kpoints)'



    # if spin[1] == "up" && spin[2] == "up"
    #     n_kpoints, n_bnds, energies_1, data_1 = load_filtered_data(FILE_PATH_BS*"_up.csv",FILE_PATH_BS*"_energies_up.csv", a_orbital[1], n_atom[1])
    #     tmp = sum(eachcol(data_1))
    #     tmp_reshape = reshape(tmp, length(n_atom[1]),:)
    #     tmp = vec(sum(tmp_reshape, dims=1))
    #     bands_1 = reshape(tmp, n_bnds, n_kpoints)'
    #     # Second bandstructure
    #     n_kpoints, n_bnds, energies_2, data_2 = load_filtered_data(FILE_PATH_BS*"_up.csv",FILE_PATH_BS*"_energies_up.csv", a_orbital[2], n_atom[2])
    #     tmp = sum(eachcol(data_2))
    #     tmp_reshape = reshape(tmp, length(n_atom[2]),:)
    #     tmp = vec(sum(tmp_reshape, dims=1))
    #     bands_2 = reshape(tmp, n_bnds, n_kpoints)'
    # elseif spin[1] == "down" && spin[2] == "down"
    #     n_kpoints, n_bnds, energies_1, data_1 = load_filtered_data(FILE_PATH_BS*"_dn.csv",FILE_PATH_BS*"_energies_dn.csv", a_orbital[1], n_atom[1])
    #     tmp = sum(eachcol(data_1))
    #     tmp_reshape = reshape(tmp, length(n_atom[1]),:)
    #     tmp = vec(sum(tmp_reshape, dims=1))
    #     bands_1 = reshape(tmp, n_bnds, n_kpoints)'
    #     # Second bandstructure
    #     n_kpoints, n_bnds, energies_2, data_2 = load_filtered_data(FILE_PATH_BS*"_dn.csv",FILE_PATH_BS*"_energies_dn.csv", a_orbital[2], n_atom[2])
    #     tmp = sum(eachcol(data_2))
    #     tmp_reshape = reshape(tmp, length(n_atom[2]),:)
    #     tmp = vec(sum(tmp_reshape, dims=1))
    #     bands_2 = reshape(tmp, n_bnds, n_kpoints)'
    # elseif spin[1] == "up" && spin[2] == "down"
    #     n_kpoints, n_bnds, energies_1, data_1 = load_filtered_data(FILE_PATH_BS*"_up.csv",FILE_PATH_BS*"_energies_up.csv", a_orbital[1], n_atom[1])
    #     tmp = sum(eachcol(data_1))
    #     tmp_reshape = reshape(tmp, length(n_atom[1]),:)
    #     tmp = vec(sum(tmp_reshape, dims=1))
    #     bands_1 = reshape(tmp, n_bnds, n_kpoints)'
    #     # Second bandstructure
    #     n_kpoints, n_bnds, energies_2, data_2 = load_filtered_data(FILE_PATH_BS*"_dn.csv",FILE_PATH_BS*"_energies_dn.csv", a_orbital[2], n_atom[2])
    #     tmp = sum(eachcol(data_2))
    #     tmp_reshape = reshape(tmp, length(n_atom[2]),:)
    #     tmp = vec(sum(tmp_reshape, dims=1))
    #     bands_2 = reshape(tmp, n_bnds, n_kpoints)'
    # else
    #     n_kpoints, n_bnds, energies_1, data_1 = load_filtered_data(FILE_PATH_BS*"_dn.csv",FILE_PATH_BS*"_energies_dn.csv", a_orbital[1], n_atom[1])
    #     tmp = sum(eachcol(data_1))
    #     tmp_reshape = reshape(tmp, length(n_atom[1]),:)
    #     tmp = vec(sum(tmp_reshape, dims=1))
    #     bands_1 = reshape(tmp, n_bnds, n_kpoints)'
    #     # second bandstructure
    #     n_kpoints, n_bnds, energies_2, data_2 = load_filtered_data(FILE_PATH_BS*"_up.csv",FILE_PATH_BS*"_energies_up.csv", a_orbital[2], n_atom[2])
    #     tmp = sum(eachcol(data_2))
    #     tmp_reshape = reshape(tmp, length(n_atom[2]),:)
    #     tmp = vec(sum(tmp_reshape, dims=1))
    #     bands_2 = reshape(tmp, n_bnds, n_kpoints)'
    # end


    # Calculate sections and x_axis
    secs = cumsum([0;norm.(eachrow(diff(sym_points,dims=1)),2)[1:2:end]])
    x_axis = vcat(
        range.(secs[1:end-2],secs[2:end-1],length=50)...,
        range(secs[end-1],secs[end], length = 70)...
        )
    ## Plot data with dots
    plot!(p1,x_axis,reshape(energies_1, n_bnds, n_kpoints)'.-fermi_energy,
        markersize=LINE_WIDTH[1],
        seriescolor=LINE_COLOR[1],
        seriestype = :scatter,
        seriesalpha = bands_1,
        label=false,xticks=false
        )
    plot!(p1,(x_axis[1],energies_1[1]),
        markersize=LINE_WIDTH[2],
        seriescolor=LINE_COLOR[1],
        seriestype = :scatter,
        label="spin-up"
        )
    plot!(p1,x_axis,reshape(energies_2, n_bnds, n_kpoints)' .- fermi_energy,
        markersize=LINE_WIDTH[2],
        seriescolor=LINE_COLOR[2],
        seriestype = :scatter,
        seriesalpha = bands_2,
        label=false,xticks=false
        )
    plot!(p1,(x_axis[1],energies_2[1]),
        markersize=LINE_WIDTH[2],
        seriescolor=LINE_COLOR[2],
        seriestype = :scatter,
        label="spin-down"
        )
    ## Plot data with lines

    # plot!(p1,x_axis,reshape(energies_1, n_bnds, n_kpoints)',
    #     lw=LINE_WIDTH[1],
    #     lc=LINE_COLOR[1],
    #     ls = :dot,
    #     seriesalpha = bands_1,
    #     label=false,xticks=false
    #     )

    # plot!(p1,x_axis,reshape(energies_2, n_bnds, n_kpoints)',
    #     lw=LINE_WIDTH[2],
    #     seriescolor=LINE_COLOR[2],
    #     ls = :dot,
    #     seriesalpha = bands_2,
    #     label=false,xticks=false
    #     )

    
    # # Plot fermi energy and sections
    # hline!(p1,[fermi_energy],
    #     lw=LINE_WIDTH[3],
    #     lc=LINE_COLOR[3],
    #     ls=:dashdot,label=false
    #     )
    hline!(p1,[0],
        lw=LINE_WIDTH[3],
        lc=LINE_COLOR[3],
        ls=:dashdot,label=false
    )
    vline!(p1,secs,lw=LINE_WIDTH[4],
        lc=:black,
        ls=:solid,label=false
        )

    # Set plot properties
    plot!(p1,grid=false, ylabel=L"Energy - $E_F$  (eV)",
        ylims=(-1,1),
        xlims=(0,x_axis[end]), # Change to x_axis[end]
    )
    plot!(p1, xticks = (secs,label_sym_points(FILE_PATH_SYM)),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12)
    plot!(p1, label = :topright)

    annotate!(p1,0.75,fermi_energy+ANNOTATE_POSITION, Plots.text("Fermi energy", ANNOTATE_SIZE))

    return p1
    #return x_axis, reshape(energies_1.px, n_bnds, n_kpoints), bands_1
end

# Code to plot the Projected Bandstructure fo a non-spin polarized material. The Bandstructure lines are colored (opacity)
# according to the density of states projected on the chosen orbital of the chosen atom.
# The input files are the paths to the PROCAR file, the KPOINTS file (FILE_PATH_SYM). The next inputs are a vector with the desired atoms to be plotted,
# and a vector with the desired orbitals to be plotted.
# For instance, if I want to plot the bandstructure of the d orbitals of the first and second atoms and 
# s-p orbitals of the second atom, the inputs are n_atom = [[1,2],[2]], and a_orbital = [["d"],["s","p"]].
# Note that the number labels of the n_atom are according with the POSCAR file, i.e., the first half of atoms are Nb, and the second half are P.
# The output are two plots with the weighted Projected Bandstructure of the material. The plot includes the Fermi energy and the symmetry points.
################### Parameters to tune the plot ###################
# LINE_WIDTH: Array with the width of the lines of the plot. The first element is the width of the bandstructure lines, 
# the second element is the width of the Fermi energy line, and the third element is the vertical symmetry lines.
# LINE_COLOR: Array with the color of the lines of the plot. The first element is the color of the spin-up bandstructure lines,
# the second element is the color of the spin-down bandstructure lines, and the third element is the color of Fermi line.
function Bandstructure_nospin_plot(PATH_DF::String,FILE_PATH_SYM::String,n_atom::Vector{Vector{Int64}},a_orbital::Vector{Vector{String}})
    p1 = plot();
        # Parameters for the plot
    # LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:red, :blue, :black]
    ANNOTATE_POSITION = 0.1
    ANNOTATE_SIZE = 8

    fermi_energy = CSV.File(open(PATH_DF*"df_energies_total.csv"); limit=5).energy[5]
    sym_points = readdlm(FILE_PATH_SYM)[5:end,1:3]
    
    n_kpoints, n_bnds, energies_1, data_1 = load_filtered_data(PATH_DF*"df_bands_total.csv", PATH_DF*"df_energies_total.csv", a_orbital[1], n_atom[1])
    tmp = sum(eachcol(data_1))
    tmp_reshape = reshape(tmp, length(n_atom[1]),:)
    tmp = vec(sum(tmp_reshape, dims=1))
    bands_1 = reshape(tmp, n_bnds, n_kpoints)'
    # Second bandstructure
    n_kpoints, n_bnds, energies_2, data_2 = load_filtered_data(PATH_DF*"df_bands_total.csv", PATH_DF*"df_energies_total.csv", a_orbital[2], n_atom[2])
    tmp = sum(eachcol(data_2))
    tmp_reshape = reshape(tmp, length(n_atom[2]),:)
    tmp = vec(sum(tmp_reshape, dims=1))
    bands_2 = reshape(tmp, n_bnds, n_kpoints)'

    # Calculate sections and x_axis
    secs = cumsum([0;norm.(eachrow(diff(sym_points,dims=1)),2)[1:2:end]])
    x_axis = vcat(range.(secs[1:end-2],secs[2:end-1],length=50)...,range(secs[end-1],secs[end], length = 70)...)

    # Plot data

    # plot!(p1,x_axis,reshape(energies_1, n_bnds, n_kpoints)',
    #     lw=LINE_WIDTH[1],
    #     lc=LINE_COLOR[1],
    #     linealpha=bands_1,
    #     ls=:dot,
    #     label=false,xticks=false
    #     )


    # plot!(p1,x_axis,reshape(energies_2, n_bnds, n_kpoints)',
        # lw=LINE_WIDTH[2],
        # lc=LINE_COLOR[2],
        # linealpha=bands_2,
        # ls=:dot,
        # label=false,xticks=false
        # )


    plot!(p1,x_axis,reshape(energies_1, n_bnds, n_kpoints)'.-fermi_energy,
        markersize=LINE_WIDTH[1],
        seriescolor=LINE_COLOR[1],
        seriestype = :scatter,
        seriesalpha = bands_1,
        label=false,xticks=false
        )
    plot!(p1,(x_axis[1],energies_1[1]),
        markersize=LINE_WIDTH[2],
        seriescolor=LINE_COLOR[1],
        seriestype = :scatter,
        label="d orbitals"
        )

    plot!(p1,x_axis,reshape(energies_2, n_bnds, n_kpoints)'.-fermi_energy,
        markersize=LINE_WIDTH[2],
        seriescolor=LINE_COLOR[2],
        seriestype = :scatter,
        seriesalpha = 20*bands_2, ############### Change the value of the opacity
        label=false,xticks=false
        )
    plot!(p1,(x_axis[1],energies_2[1]),
        markersize=LINE_WIDTH[2],
        seriescolor=LINE_COLOR[2],
        seriestype = :scatter,
        label="s-p orbitals"
        )
    
    # Plot fermi energy and sections
    hline!(p1,[0],
        lw=LINE_WIDTH[3],
        lc=LINE_COLOR[3],
        ls=:dashdot,label=false
        )
    vline!(p1,secs,lw=LINE_WIDTH[4],
        lc=:black,
        ls=:solid,label=false
        )

    # Set plot properties
    plot!(p1,grid=false, ylabel=L"Energy - $E_F$  (eV)",
        ylims=(-3,3),
        xlims=(0,x_axis[end]),
    )
    plot!(p1, xticks = (secs,label_sym_points(FILE_PATH_SYM)),
        ytickfontsize  = 10,
        xtickfontsize  = 10,
        yguidefontsize = 18,
        legendfontsize = 12
        )
    annotate!(p1,1.0,fermi_energy+ANNOTATE_POSITION, Plots.text("Fermi energy", ANNOTATE_SIZE))

    #return x_axis, reshape(energies_1.px, n_bnds, n_kpoints), bands_1
    return p1
end


function Total_Bandstructure_plot(FILE_PATH_BS::String,FILE_PATH_SYM::String)
    p1 = plot();
     # Parameters for the plot
    LINE_WIDTH = [1.75, 0.75, 1.0]
    LINE_COLOR = [:red, :blue, :black]
    ANNOTATE_POSITION = 0.2
    ANNOTATE_SIZE = 8
 
    fermi_energy = CSV.File(open(FILE_PATH_BS*"_energies_up.csv"); limit=5).energy[5]
    sym_points = readdlm(FILE_PATH_SYM)[5:end,1:3]

    #Calculation of the energies
    n_kpoints1, n_bnds1, energies_1, _ = load_filtered_data(FILE_PATH_BS*"_up.csv", FILE_PATH_BS*"_energies_up.csv", ["s"], [1])
    n_kpoints2, n_bnds2, energies_2, _ = load_filtered_data(FILE_PATH_BS*"_dn.csv", FILE_PATH_BS*"_energies_dn.csv", ["s"], [2])


    # Calculate sections and x_axis
    secs = cumsum([0;norm.(eachrow(diff(sym_points,dims=1)),2)[1:2:end]])
    x_axis = vcat(range.(secs[1:end-1],secs[2:end],length=20)...)
    # Plot data
    plot!(p1,x_axis,reshape(energies_1, n_bnds1, n_kpoints1)',
        lw=LINE_WIDTH[1],
        lc=LINE_COLOR[1],
        linealpha=0.8,
        ls=:dot,
        label=false,xticks=false
        )

    plot!(p1,x_axis,reshape(energies_2, n_bnds2, n_kpoints2)',
        lw=LINE_WIDTH[1],
        lc=LINE_COLOR[2],
        linealpha= 0.8,
        ls=:dot,
        label=false,xticks=false
        )

    
    # Plot fermi energy and sections
    hline!(p1,[fermi_energy],
        lw=LINE_WIDTH[2],
        lc=LINE_COLOR[3],
        ls=:dashdot,label=false
        )
    vline!(p1,secs,lw=LINE_WIDTH[3],
        lc=:black,
        ls=:solid,label=false
        )

    # Set plot properties
    # plot!(p1,grid=false, ylabel="Energy (eV)", ylim=(-2.5,3.5))
    plot!(p1,grid=false, ylabel="Energy (eV)",
        ylims=(fermi_energy-1,fermi_energy+1),
        xlims=(0,x_axis[end]),
    )
    plot!(p1, xticks = (secs,label_sym_points(FILE_PATH_SYM)))
    annotate!(p1,1.5,fermi_energy+ANNOTATE_POSITION, Plots.text("Fermi energy", ANNOTATE_SIZE))
    display(p1)
end


function Total_nBandstructure_plot(FILE_PATH_BS::String,FILE_PATH_SYM::String)
    p1 = plot();
     # Parameters for the plot
    LINE_WIDTH = [1.75, 0.75, 1.0]
    LINE_COLOR = [:red, :blue, :black]
    ANNOTATE_POSITION = 0.2
    ANNOTATE_SIZE = 8
 
    fermi_energy = CSV.File(open(FILE_PATH_BS*"df_energies_total.csv"); limit=5).energy[5]
    tot_atoms = CSV.File(open(FILE_PATH_BS*"df_energies_total.csv"); limit=4).energy[4]

    sym_points = readdlm(FILE_PATH_SYM)[5:end,1:3]

    #Calculation of the energies
    n_kpoints, n_bnds, energies_2, data_2 = load_filtered_data(FILE_PATH_BS*"df_bands_total.csv", FILE_PATH_BS*"df_energies_total.csv", ["s","p","d"], [0])
    tmp = sum(eachcol(data_2))
    tmp_reshape = reshape(tmp, length(tot_atoms),:)
    tmp = vec(sum(tmp_reshape, dims=1))
    bands_2 = reshape(tmp, n_bnds, n_kpoints)'

    # Calculate sections and x_axis
    secs = cumsum([0;norm.(eachrow(diff(sym_points,dims=1)),2)[1:2:end]])
    x_axis = vcat(range.(secs[1:end-2],secs[2:end-1],length=50)...,range(secs[end-1],secs[end], length = 70)...)
    # Plot data
    plot!(p1,x_axis,reshape(energies_2, n_bnds, n_kpoints)'.-fermi_energy,
        markersize=LINE_WIDTH[2],
        seriescolor=LINE_COLOR[2],
        seriestype = :scatter,
        seriesalpha = bands_2,
        label=false,xticks=false
        )
    plot!(p1,(x_axis[1],energies_2[1]),
        markersize=LINE_WIDTH[2],
        seriescolor=LINE_COLOR[2],
        seriestype = :scatter,
        label="VASP"
        )
    
    # Plot fermi energy and sections
    hline!(p1,[fermi_energy],
        lw=LINE_WIDTH[2],
        lc=LINE_COLOR[3],
        ls=:dashdot,label=false
        )
    vline!(p1,secs,lw=LINE_WIDTH[3],
        lc=:black,
        ls=:solid,label=false
        )

    # Set plot properties
    # plot!(p1,grid=false, ylabel="Energy (eV)", ylim=(-2.5,3.5))
    plot!(p1,grid=false, ylabel="Energy (eV)",
        ylims=(-1,1),
        xlims=(0,x_axis[end]),
    )
    plot!(p1, xticks = (secs,label_sym_points(FILE_PATH_SYM)))
    annotate!(p1,1.5,fermi_energy+ANNOTATE_POSITION, Plots.text("Fermi energy", ANNOTATE_SIZE))
    display(p1)
end


# Surface plot of the bandstructure
# Bandstructure_spin_plot(FILE_BAND1,FILE_SYM1,[[3,6,9,12],[3,6,9,12]],["up","down"],[["d"],["d"]])
# Bandstructure_spin_plot(FILE_BAND1,FILE_SYM1,[[51,54,57,60],[51,54,57,60]],["up","down"],[["s","p"],["s","p"]])

# Bandstructure_spin_plot(FILE_BAND1,FILE_SYM1,[[13,16,19,22,25,28,31,34,37,40,43,46],[13,16,19,22,25,28,31,34,37,40,43,46]],["up","down"],[["d"],["d"]])
# Bandstructure_spin_plot(FILE_BAND1,FILE_SYM1,[collect(49:3:96),collect(49:3:96)],["up","down"],[["s","p"],["s","p"]])


# Bandstructure_nospin_plot(FILE_BAND,FILE_SYM,[[3,6,9,12],[51,54,57,60]],[["d"],["s","p"]])
# Bandstructure_nospin_plot(FILE_BAND,FILE_SYM,[[13,16,19,22,25,28,31,34,37,40,43,46],collect(49:3:96)],[["d"],["s","p"]])
