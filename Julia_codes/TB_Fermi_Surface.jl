using DelimitedFiles
using CSV
using DataFrames
using LinearAlgebra
using FileIO
using SparseArrays
using Arpack
using HDF5
using JLD2

const PATH_DF = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/df/"

function closest_elements(M::Array{T}, m::T, n::Int=10) where T
    # Calculate the absolute differences
    differences = abs.(M .- m)
    
    # Get the indices of the sorted differences
    sorted_indices = sortperm(differences)
    
    # Select the first n elements based on the sorted indices
    closest_elements = M[sorted_indices[1:n]]
    
    return closest_elements .- m, sorted_indices[1:n]
end

# Function to compute the tight-binding Hamiltonian
function tb_Ham(PATH_DF::String, kx::Float64, ky::Float64)
    alat = 6.3427 # Lattice parameter in Å
    num_wan = h5open(joinpath(PATH_DF, "NbP_TB.h5"), "r") do file
        size(read(file,"unitcell_0_0_0/Hamiltonian"))[1] # Load the matrix data
    end
    H_k =  zeros(Complex{Float64}, num_wan, num_wan) # Initialize Hamiltonian matrix
    # Loop over the .jld2 files
    h5open(joinpath(PATH_DF, "NbP_TB.h5"), "r") do file
        # Iterate over each group in the file
        for group_name in keys(file)
            Rx, Ry, Rz = parse.(Int, split(group_name, r"[_\.]")[2:4])
            if Rz == 0
                # Update the Hamiltonian matrix with the corresponding exponential term
                hamiltonian = read(file,"$(group_name)/Hamiltonian") # Load the matrix data
                H_k += hamiltonian * exp(im*(Rx*kx + Ry*ky))
            end

            # Update the Hamiltonian matrix with the corresponding exponential term
            # hamiltonian = read(file,"$(group_name)/Hamiltonian") # Load the matrix data
            # H_k += hamiltonian * exp(im*(Rx*kx + Ry*ky))

        end
    end
    # Return the Hamiltonian matrix
    return H_k # Returns the nearest neighbor interactions from the d orbitals to the p orbitals
end



# Function to project the plane
function plane_projector(PATH_DF::String,Nx::Int64,Ny::Int64, e_fermi::Float64)
    a_l = 3.335 # 2D lattice square parameter a_l = b_l = 3.335 Å = 3.335*10^(-10) m
    proj_data_dict = Dict{Tuple{Int64,Int64}, Dict{String, Any}}()
    # proj_data_dict = Dict{Tuple{Int64, Int64}, Dict{String, Any}}()
    kx = range(-pi, pi, length = Nx) # Range of ky values
    ky = range(-pi, pi, length = Ny) # Range of kz values
    k = 10;

    eigenvalues_array = Array{Float64, 3}(undef, Nx, Ny, k)
    two_lay_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    one_lay_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    second_lay_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    # d- orbital projections
    d1_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    d2_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    d3_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    d4_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    d5_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    # sp- orbital projections
    sp1_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    sp2_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    sp3_proj = Array{Float64, 3}(undef, Nx, Ny, k)
    sp4_proj = Array{Float64, 3}(undef, Nx, Ny, k)


    # Threads.@threads for ix in 1:Nx
    for ix in 1:Nx
	    for iy in 1:Ny
		    # Compute the Hamiltonian for the given kx, ky pair
		    H = tb_Ham(PATH_DF, kx[ix], ky[iy])

		    # # Perform eigenvalue decomposition
		    # eigvals, eigenvecs = eigs(H, nev = k, which=:LM, sigma= e_fermi)
		    eigvals, eigenvecs = eigen(H)
            # p = sortperm(real(eigvals))
            # eigvals = eigvals[p]
            eigenvalues_array[ix, iy, :], sorted_indices = closest_elements(real(eigvals),  e_fermi)
            eigenvecs = eigenvecs[:, sorted_indices]
            # eigenvecs = eigenvecs[:, p]
		    # Calculate the sum of squared magnitudes for the first 36 elements of each eigenvector
            one_lay_proj[ix, iy, :] = [sum(abs.(eigenvecs[1:20, i]).^2) for i in 1:k]
            second_lay_proj[ix, iy, :] = [sum(abs.(eigenvecs[21:36, i]).^2) for i in 1:k]
            two_lay_proj[ix, iy, :] = [sum(abs.(eigenvecs[1:36, i]).^2) for i in 1:k]
            # Calculate the sum of squared magnitudes of each surface d orbital
            d1_proj[ix, iy, :] = [sum(abs.(eigenvecs[1:5:20, i]).^2) for i in 1:k]
            d2_proj[ix, iy, :] = [sum(abs.(eigenvecs[2:5:20, i]).^2) for i in 1:k]
            d3_proj[ix, iy, :] = [sum(abs.(eigenvecs[3:5:20, i]).^2) for i in 1:k]
            d4_proj[ix, iy, :] = [sum(abs.(eigenvecs[4:5:20, i]).^2) for i in 1:k]
            d5_proj[ix, iy, :] = [sum(abs.(eigenvecs[5:5:20, i]).^2) for i in 1:k]
            # Calculate the sum of squared magnitudes of each surface sp orbital
            sp1_proj[ix, iy, :] = [sum(abs.(eigenvecs[21:4:36, i]).^2) for i in 1:k]
            sp2_proj[ix, iy, :] = [sum(abs.(eigenvecs[22:4:36, i]).^2) for i in 1:k]
            sp3_proj[ix, iy, :] = [sum(abs.(eigenvecs[23:4:36, i]).^2) for i in 1:k]
            sp4_proj[ix, iy, :] = [sum(abs.(eigenvecs[24:4:36, i]).^2) for i in 1:k]

            # Store the projections in the dictionary
            proj_data_dict[(ix, iy)] = Dict(
                "2-layers"     => two_lay_proj[ix, iy, :],
                "1-layer"      => one_lay_proj[ix, iy, :],
                "second-layer" => second_lay_proj[ix, iy, :],
                "d1"  => d1_proj[ix, iy, :],
                "d2"  => d2_proj[ix, iy, :],
                "d3"  => d3_proj[ix, iy, :],
                "d4"  => d4_proj[ix, iy, :],
                "d5"  => d5_proj[ix, iy, :],
                "sp1" => sp1_proj[ix, iy, :],
                "sp2" => sp2_proj[ix, iy, :],
                "sp3" => sp3_proj[ix, iy, :],
                "sp4" => sp4_proj[ix, iy, :]
            )

	    end
    end

    # Return the projected plane
    return eigenvalues_array, proj_data_dict
end

# Usage
# save_data_to_hdf5(eigvalues_array, eigen_data_dict, "eigen_data.h5")
############################################################################################################
# eigvalues_array, eigen_data_dict = plane_projector(@__DIR__, 20, 20, 10.3786)
# eigvalues_array, eigen_data_dict = plane_projector(PATH_TB, 200, 200, 10.3786);

# @save "eigenvalues_array.jld2" eigvalues_array
# @save "eigenvectors_dict.jld2" eigen_data_dict

