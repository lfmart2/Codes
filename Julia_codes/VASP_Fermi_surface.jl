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

# The following code reads the PROCAR file, the DOSCAR file, and the KPOINTS file.
const FILE_EIGENVAL = raw"/home/r_floren/BandStructure/VASP/NbP/Supercell_2_2_3/Spin-QE/DOS.spin.QE/EIGENVAL"


function filter_file(FILE_PATH_BAND::String, FILE_PATH_FE::String)
    df = open(FILE_PATH_BAND) do io
        dropmissing(CSV.read(io, DataFrame; header=["ion", "s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2", "tot"], skipto =2, ignoreemptyrows=true, delim=' ', ignorerepeated=true, footerskip=1), :ion)
    end
    df = deleteat!(df,findall(i -> (i == "k-point") || (i == "ion"), df.ion))
    if length(findall(i->i=="#", df.ion)) == 1
        init_info = findall(i->i=="#", df.ion)[1]
        info = DataFrame(band = ["Non-polarized", "#-kpoints:","#-bands:","#-ions:","E-Fermi:"],
            energy = ["",parse(Int64,df[init_info,4]),parse(Int64,df[init_info,8]),df[init_info,12],extract_fermi_energy(FILE_PATH_FE)])
        energies = parse.(Float64,df.px[findall(i -> i=="band", df.ion)])
        tmp = DataFrame("band" => repeat(1:info.energy[3],info.energy[2]), "energy" => energies)
        df = select(deleteat!(df,findall(i -> i=="#" || i == "band", df.ion)),collect(1:11))
        df.ion[findall(i -> i=="tot", df.ion)] .= "0"
        col = names(df)
        df[!,col[1]] .= parse.(Int64,df[!,col[1]])
        for i in 2:length(col)
            df[!,col[i]] .= parse.(Float64,df[!,col[i]])
        end
        # .csv file for writing with name POSCAR.csv
        CSV.write(FILE_PATH_BAND*"_bnd.csv", df)
        
        CSV.write(FILE_PATH_BAND*"_energies.csv", vcat(info,tmp))
    else
        spn_up_info = findall(i->i=="#", df.ion)[1]
        spn_dn_info = findall(i->i=="#", df.ion)[2]
        info_up = DataFrame(band = ["S-up:", "#-kpoints:","#-bands:","#-ions:","E-Fermi:"],
            energy = ["",parse(Int64,df[spn_up_info,4]),parse(Int64,df[spn_up_info,8]),df[spn_up_info,12],extract_fermi_energy(FILE_PATH_FE)])
        info_dn = DataFrame(band = ["S-up:", "#-kpoints:","#-bands:","#-ions:","E-Fermi:"],
            energy = ["",parse(Int64,df[spn_dn_info,4]),parse(Int64,df[spn_dn_info,8]),df[spn_dn_info,12],extract_fermi_energy(FILE_PATH_FE)])
        energies_up = parse.(Float64,df.px[findall(i -> i=="band", df.ion)])[1:end÷2]
        energies_dn = parse.(Float64,df.px[findall(i -> i=="band", df.ion)])[end÷2+1:end]
        tmp_up = DataFrame("band" => repeat(1:info_up.energy[3],info_up.energy[2]), "energy" => energies_up)
        tmp_dn = DataFrame("band" => repeat(1:info_dn.energy[3],info_dn.energy[2]), "energy" => energies_dn)
        
        df = select(deleteat!(df,findall(i -> i=="#" || i == "band", df.ion)),collect(1:11))
        df.ion[findall(i -> i=="tot", df.ion)] .= "0"
        col = names(df)
        df[!,col[1]] .= parse.(Int64,df[!,col[1]])
        for i in 2:length(col)
            df[!,col[i]] .= parse.(Float64,df[!,col[i]])
        end
        # .csv file for writing with name POSCAR.csv
        CSV.write(FILE_PATH_BAND*"_up.csv", df[1:end÷2,:])
        CSV.write(FILE_PATH_BAND*"_dn.csv", df[end÷2+1:end,:])

        CSV.write(FILE_PATH_BAND*"_energies_up.csv", vcat(info_up,tmp_up))
        CSV.write(FILE_PATH_BAND*"_energies_dn.csv", vcat(info_dn,tmp_dn))
    end
    
end


function Fermi_surface(FILE_PATH_EIG::String)
    p1 = plot();
        # Parameters for the plot
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:red, :blue, :black]
    ANNOTATE_POSITION = 0.05
    ANNOTATE_SIZE = 9

    fermi_energy = 7.24621144
    num_bands = parse(Int64,split(read_line(FILE_PATH_EIG,6))[end])


    df = open(FILE_PATH_EIG) do io
        dropmissing(CSV.read(io, DataFrame; header=["band", "En_up", "En_dn"], delim=' ',ignorerepeated=true, skipto = 8)[:,1:3],:band)
    end

    # Chunk the data into spin-up and spin-down
    k_vecs = df[1:num_bands+1:end,:]
    tmp0 = 1:552+1:length(df.band)
    println(tmp0)
    tmp = zeros(4,4,num_bands)
    tmp1 = 1
    for i in 1:4
        for j in 1:4
            tmp[i,j,:] = df.En_up[tmp0[tmp1]+1:tmp0[tmp1]+num_bands]
            tmp1 += 1
        end
    end
    return tmp
    #return x_axis, reshape(energies_1.px, n_bnds, n_kpoints), bands_1
end
# p1 = plot()

# contour!(p1,kx,ky,tmp1[:,:,435].-7.24621144)

# contour!(p1,-3:1:0,ky,reverse(tmp1[:,:,435], dims=2).-7.24621144)

# contour!(p1,kx,-3:1:0,reverse(tmp1[:,:,435], dims=1).-7.24621144)

# contour!(p1,-3:1:0,-3:1:0,reverse(reverse(tmp1[:,:,435], dims=1),dims=2).-7.24621144)