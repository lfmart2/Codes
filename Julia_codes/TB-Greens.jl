using Plots
using DelimitedFiles
using HDF5
using CSV
using Plots.PlotMeasures
using DataFrames
using LinearAlgebra
using LaTeXStrings

# Define the path to the directory containing the DataFrames
const PATH_TB = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/df/"
const PATH_H5 = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/df/Planes/"



function load_planes(DF_PATH::String)
    h5open(DF_PATH*"Planes/NbP_layers.h5", "w") do file
        for l in 1:12
            group_name = "layer_$(l)"
            g = create_group(file, group_name)
            for csv_file in readdir(DF_PATH)[endswith.(readdir(DF_PATH),".csv")]
                df = CSV.read(joinpath(DF_PATH, csv_file), DataFrame)
                layer_size = 36 # Number of Wannier functions
                layer = zeros(Complex{Float64}, layer_size, layer_size) # Initialize Hamiltonian matrix
                key = [df[1,:].Rx, df[1,:].Ry, df[1,:].Rz]
                for i in 1:layer_size
                    for j in 1:layer_size
                        tmp = filter([:n₀, :nᵢ] => (n₀, nᵢ) -> n₀ == 36*(l-1)+i && nᵢ == 36*(l-1)+j, df)
                        layer[i, j] += Complex(tmp.En[1], tmp.iEn[1])
                    end
                end
                g["pos_$(key[1])_$(key[2])_$(key[3])"] = layer  # Directly save the dense matrix as a dataset
            end
        end
    end
end

function load_planes_int(DF_PATH::String)
    h5open(DF_PATH*"Planes/NbP_layers_int.h5", "w") do file
        for l in 1:11
            for ll in l+1:12
                group_name = "interaction_$(l)_to_$(ll)"
                g = create_group(file, group_name)
                for csv_file in readdir(DF_PATH)[endswith.(readdir(DF_PATH),".csv")]
                    df = CSV.read(joinpath(DF_PATH, csv_file), DataFrame)
                    layer_size = 36 # Number of Wannier functions
                    layer = zeros(Complex{Float64}, layer_size, layer_size) # Initialize Hamiltonian matrix
                    key = [df[1,:].Rx, df[1,:].Ry, df[1,:].Rz]
                    for i in 1:layer_size
                        for j in 1:layer_size
                            tmp = filter([:n₀, :nᵢ] => (n₀, nᵢ) -> n₀ == 36*(l-1)+i && nᵢ == 36*(ll-1)+j, df)
                            layer[i, j] += Complex(tmp.En[1], tmp.iEn[1])
                        end
                    end
                    g["pos_$(key[1])_$(key[2])_$(key[3])"] = layer  # Directly save the dense matrix as a dataset
                end
            end
        end
    end
end

function H_k(PATH_H5::String, kx::Float64, ky::Float64, lay::Int)
    num_wan = h5open(joinpath(PATH_H5,"NbP_layers.h5"), "r") do file
        size(read(file,"layer_1/pos_0_0_0"))[1] # Load the matrix data
    end
    Hₖ =  zeros(Complex{Float64}, num_wan, num_wan) # Initialize Hamiltonian matrix
    # Loop over the hdf5 file and load the Hamiltonian matrices
    for file in h5read(joinpath(PATH_H5,"NbP_layers.h5"), "layer_$(lay)")
        # Iterate over each group in the file
        Rx, Ry, Rz = parse.(Int, split(file[1], r"[_\.]")[2:4])
        if Rz == 0
            # Update the Hamiltonian matrix with the corresponding exponential term
            # hamiltonian = file[2] # Load the matrix data
            Hₖ += file[2] * exp(im*(Rx*kx + Ry*ky))
        end
    end
    return Hₖ
end


function T_k(PATH_H5::String, kx::Float64, ky::Float64, lay₀::Int, layᵢ::Int)
    num_wan = h5open(joinpath(PATH_H5,"NbP_layers_int.h5"), "r") do file
        size(read(file,"interaction_1_to_2/pos_0_0_0"))[1] # Load the matrix data
    end
    Hₖ =  zeros(Complex{Float64}, num_wan, num_wan) # Initialize Hamiltonian matrix
    # Loop over the hdf5 file and load the Hamiltonian matrices
    for file in h5read(joinpath(PATH_H5,"NbP_layers_int.h5"), "interaction_$(lay₀)_to_$(layᵢ)")
        # Iterate over each group in the file
        Rx, Ry, Rz = parse.(Int, split(file[1], r"[_\.]")[2:4])
        if Rz == 0
            # Update the Hamiltonian matrix with the corresponding exponential term
            Hₖ += file[2] * exp(im*((Rx + 0.25*isodd(lay₀))*kx + (0.25*iseven(lay₀) + Ry)*ky))
        end
    end
    return Hₖ # Returns the neares neighbor interactions from the d orbitals to the p orbitals
end



function bulk_Greens(PATH_H5::String, kx::Float64, ky::Float64, ω::ComplexF64)
    H0_bulk = H_k(PATH_H5, kx, ky, 12)
    T_w = T_k(PATH_H5, kx, ky,10,11)
    T_v = T_k(PATH_H5, kx, ky,11,12) 
    num_wan = size(H0_bulk)[1]
    tmp = ω*I(num_wan) .- H0_bulk .+ T_v'*inv(ω*I(num_wan) - H0_bulk)T_v .+ T_w'*inv(ω*I(num_wan) - H0_bulk)T_w
    return (tmp .+ sqrt( tmp^2 - 4*T_w'*T_w ) ) ./2

end

function surf_Greens(PATH_H5::String, kx::Float64, ky::Float64, ω::ComplexF64)
    num_wan = size(H_k(PATH_H5, 0.0, 0.0, 12))[1]
    tmp     = ω *I(num_wan) .- H_k(PATH_H5,kx,ky,9) .- T_k(PATH_H5, kx, ky,10,11)' * bulk_Greens(PATH_H5, kx, ky, ω) * T_k(PATH_H5, kx, ky,10,11)
    for i in 9:-1:1
        tmp = ω *I(num_wan) .- H_k(PATH_H5,kx,ky,i) .- T_k(PATH_H5, kx, ky,i,i+1)' * tmp * T_k(PATH_H5, kx, ky,i,i+1)
    end
    return -imag(inv(tmp))/pi
end

function plot_FS(PATH_H5::String, Nx::Int64, Ny::Int64, e_fermi::Float64)
    kx = range(-π, stop=π, length=Nx)
    ky = range(-π, stop=π, length=Ny)
    FS = zeros(Float64, Nx, Ny)
    for i in 1:Nx
        for j in 1:Ny
            FS[i,j] = sum(diag(surf_Greens(PATH_H5, kx[i], ky[j], Complex(e_fermi, 0.0000001))))
        end
    end
    return FS #heatmap(kx, ky, FS, c=:viridis, xlabel=L"k_x", ylabel=L"k_y", title="Fermi Surface")
end