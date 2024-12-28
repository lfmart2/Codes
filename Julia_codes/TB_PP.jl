using CSV
using Plots
using DataFrames
using LinearAlgebra
using LaTeXStrings
using FileIO
using JLD2
using DelimitedFiles
using SparseArrays
using Arpack
using HDF5


# Define the path to the directory containing the DataFrames
const PATH_TB = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/Fermi_Surface/Fermi_Surface_jld/"
# const PATH_DF = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/df/"



function plot_Fermi_surface(PATH_JLD2::String)
    p1 = plot()
    Energies = load(joinpath(PATH_JLD2, "eigenvalues_array.jld2"),"eigvalues_array")
    proj_data_dict = load(joinpath(PATH_JLD2, "projections_dict.jld2"),"proj_data_dict")
    # proj = reduce(hcat,hcat(load(tmp,"proj")...))

    length_kx = size(Energies)[1]
    length_ky = size(Energies)[2]
    k = size(Energies)[3]

    # Initialize arrays to store reconstructed data
    reconstructed_one_lay_proj    = zeros(length_kx, length_ky, k)
    reconstructed_second_lay_proj = zeros(length_kx, length_ky, k)
    reconstructed_two_lay_proj    = zeros(length_kx, length_ky, k)
    reconstructed_d1_proj  = zeros(length_kx, length_ky, k)
    reconstructed_d2_proj  = zeros(length_kx, length_ky, k)
    reconstructed_d3_proj  = zeros(length_kx, length_ky, k)
    reconstructed_d4_proj  = zeros(length_kx, length_ky, k)
    reconstructed_d5_proj  = zeros(length_kx, length_ky, k)
    reconstructed_sp1_proj = zeros(length_kx, length_ky, k)
    reconstructed_sp2_proj = zeros(length_kx, length_ky, k)
    reconstructed_sp3_proj = zeros(length_kx, length_ky, k)
    reconstructed_sp4_proj = zeros(length_kx, length_ky, k)

    # Reconstruct the arrays by accessing data from proj_data_dict
    for ix in 1:length_kx
        for iy in 1:length_ky
            proj_dict = proj_data_dict[(ix, iy)]
            
            # Fill in each array with data from the dictionary
            reconstructed_two_lay_proj[ix, iy, :] = proj_dict["2-layers"]
            reconstructed_one_lay_proj[ix, iy, :] = proj_dict["1-layer"]
            reconstructed_second_lay_proj[ix, iy, :] = proj_dict["second-layer"]
            reconstructed_d1_proj[ix, iy, :] = proj_dict["d1"]
            reconstructed_d2_proj[ix, iy, :] = proj_dict["d2"]
            reconstructed_d3_proj[ix, iy, :] = proj_dict["d3"]
            reconstructed_d4_proj[ix, iy, :] = proj_dict["d4"]
            reconstructed_d5_proj[ix, iy, :] = proj_dict["d5"]
            reconstructed_sp1_proj[ix, iy, :] = proj_dict["sp1"]
            reconstructed_sp2_proj[ix, iy, :] = proj_dict["sp2"]
            reconstructed_sp3_proj[ix, iy, :] = proj_dict["sp3"]
            reconstructed_sp4_proj[ix, iy, :] = proj_dict["sp4"]
        end
    end

    for j in 1:k
    # j=193
        contour!(p1, range(-pi, pi, length = length_ky), range(-pi, pi, length = length_ky), Energies[:,:,j],
            levels=[-0.0713999], color = 1,
            seriesalpha = reconstructed_two_lay_proj[:,:,j],
            aspect_ratio=:equal,cbar=false
            )
        contour!(p1, range(pi, 3*pi, length = length_ky), range(-pi, pi, length = length_ky), Energies[:,:,j],
            levels=[-0.0713999],color=1,
            aspect_ratio=:equal,cbar=false
            )
    end
    plot!(p1,xlabel = "kx/a", ylabel = "ky/a")
    return p1
end

function tb_proj_bandstructure(PATH_DF::String)
    e_fermi = 10.3786 # Fermi energy
    kx = vcat(zeros(50),range(0, pi, length = 50),zeros(50) .+ pi,range(pi, 0, length = 70)) # Range of kx values
    ky = vcat(range(pi, 0, length = 50), zeros(50),range(0, pi, length = 50),range(pi, 0, length = 70)) # Range of ky values
    k = 40
    # Compute the Hamiltonian for the given kx, ky pair
    H = tb_Ham1(PATH_DF, kx[1], ky[1])

    # Perform eigenvalue decomposition
    eigvals, eigenvecs = eigs(H, nev = k, which=:LM, sigma= e_fermi)

    Energies = real(eigvals)
    Projections = [sum(abs.(eigenvecs[1:36, i]).^2) for i in 1:length(eigvals)]
    for i in 2:length(kx)
        # Compute the Hamiltonian for the given kx, ky pair
        H = tb_Ham1(PATH_DF, kx[i], ky[i])

        # Perform eigenvalue decomposition
        eigvals, eigenvecs = eigs(H, nev = k, which=:LM, sigma= e_fermi)

        # Store the eigenvalues
        Energies = hcat(Energies, real(eigvals))
        Projections = hcat(Projections, [sum(abs.(eigenvecs[1:36, i]).^2) for i in 1:length(eigvals)])

    end
    return kx, ky, Energies, Projections
end




function H_RPA(PATH_DF::String)
    h5open(joinpath(PATH_DF,"NbP_TB.h5"), "w") do file
        for csv_file in readdir(PATH_DF)[endswith.(readdir(PATH_DF),".csv")]
            df = CSV.read(joinpath(PATH_DF, csv_file), DataFrame)
            num_wan = maximum(df.n₀) # Number of Wannier functions
            hamiltonian = zeros(Complex{Float64}, num_wan, num_wan) # Initialize Hamiltonian matrix
            key = [df[1,:].Rx, df[1,:].Ry, df[1,:].Rz]
            for i in 1:num_wan^2
                hamiltonian[df[i,:].n₀, df[i,:].nᵢ] += Complex(df[i,:].En, df[i,:].iEn)
            end
            group_name = "unitcell_$(key[1])_$(key[2])_$(key[3])"
            g = create_group(file, group_name)
            g["Hamiltonian"] = hamiltonian  # Directly save the dense matrix as a dataset
        end
    end
end

# Function to project the plane
function plane_projector1(PATH_DF::String,Nx::Int64,Ny::Int64, e_fermi::Float64)
    a_l = 3.335 # 2D lattice square parameter a_l = b_l = 3.335 Å = 3.335*10^(-10) m
    eigen_data_dict = Dict{Tuple{Int, Int}, Dict{String, Any}}()
    kx = range(-pi, pi, length = Nx) # Range of ky values
    ky = range(-pi, pi, length = Ny) # Range of kz values
    k = 10

    eigenvalues_array = Array{Float64, 3}(undef, Nx, Ny, k)
    for kₓ in eachindex(kx)

	    for k_y in eachindex(ky)
		    # Compute the Hamiltonian for the given kx, ky pair
		    H = tb_Ham1(PATH_DF, kx[kₓ], ky[k_y])

		    # Perform eigenvalue decomposition
		    eigvals, eigenvecs = eigs(H, nev = k, which=:LM, sigma= e_fermi)
		    
            p = sortperm(real(eigvals))
            eigvals = eigvals[p]
            eigenvalues_array[kₓ, k_y, :] = real(eigvals) .- e_fermi
		    # Calculate the sum of squared magnitudes for the first 36 elements of each eigenvector
            # result_row_vals = [sum(abs.(F_vecs[1:36, i]).^2) for i in 1:length(F.values)]
            eigen_data_dict[(kₓ, k_y)] = Dict(
                "eigenvectors" => eigenvecs[:,p] # Store top k eigenvectors
            )

	    end
    end
    # Flatten the nested arrays and subtract the Fermi energy

    return eigenvalues_array, eigen_data_dict
end

function closest_elements(M::Array{T}, m::T, n::Int=10) where T
    # Calculate the absolute differences
    differences = abs.(M .- m)
    
    # Get the indices of the sorted differences
    sorted_indices = sortperm(differences)
    
    # Select the first n elements based on the sorted indices
    closest_elements = M[sorted_indices[1:n]]
    
    return closest_elements
end