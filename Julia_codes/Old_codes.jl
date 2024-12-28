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
# const PATH_DF = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/df/"
# Store the loaded DataFrames in a global variable
const DataFramesCache = Dict{String, DataFrame}()


function load_dataframes(DF_PATH::String)
    # Load DataFrames from the specified directory if they ar2e not already loaded
    if isempty(DataFramesCache)
        for file in readdir(DF_PATH)
		if endswith(file, ".csv")  # Adjust based on your file format
                df = CSV.read(joinpath(DF_PATH, file), DataFrame)
                #  CSV.File(joinpath(DF_PATH, file)) |> DataFrame # Adjust based on your file format.
                # For Larger Datasets: If you anticipate working with larger datasets in the future or want
                # to optimize for memory usage, CSV.File would be the better choice.
                DataFramesCache[file] = df
            end
        end
    end
end



function find_H(FILE_PATH_TB::String)
    load_dataframes(FILE_PATH_TB)
    for (filename, df) in DataFramesCache
        Vector(df[1,:]) == [0, 0, 0, 0] ? println(filename) : nothing
    end
end


function max_En(PATH_TB_DF::String)
    load_dataframes(PATH_TB_DF)
    max_En  = 0
    max_iEn = 0
    lbl_En = []
    lbl_iEn = []
    for (filename, df) in DataFramesCache
        if df.Rz[1] != 0.0
            if maximum(abs.(df.En[1:end])) > max_En
                max_En = maximum(abs.(df.En[1:end]))
                lbl_En = filename
            end
            
            if maximum(abs.(df.iEn[1:end])) > max_iEn
                max_iEn = maximum(abs.(df.iEn[1:end]))
                lbl_iEn = filename
            end
        end
    end
    return max_En, max_iEn, lbl_En, lbl_iEn
end

## function to eliminate all the z ≠ 0 values for the tight-binding Hamiltonian on the current folder
function filter_2D(PATH_DF::String)
    load_dataframes(PATH_DF)
    for (filename, df) in DataFramesCache
        if df.Rz[1] != 0
            println(filename)
            rm(PATH_DF*filename)
        end
    end
end


# Function to compute the tight-binding Hamiltonian
function tb_Ham_old(PATH_DF::String, kx::Float64, ky::Float64)
    alat = 6.3427 # Lattice parameter in Å
    load_dataframes(PATH_DF)
    num_wan = maximum(
        DataFramesCache[collect(keys(DataFramesCache))[1]].n₀
        ) # Number of Wannier functions
    H_k = zeros(Complex{Float64}, num_wan, num_wan) # Initialize Hamiltonian matrix
    # Loop over the found rows
    for (filename, df) in DataFramesCache
        # Loop over the number of Wannier functions squared
        for i in 1:num_wan^2
            # Update the Hamiltonian matrix with the corresponding exponential term
            H_k[df.n₀[i], df.nᵢ[i]] += exp(im*(df.Rx[i]*kx + df.Ry[i]*ky)) * (df.En[i] + df.iEn[i]*im)
        end
    end
    # Return the Hamiltonian matrix
    return H_k # Returns the nearest neighbor interactions from the d orbitals to the p orbitals
end



# Function to project the plane
function plane_projector_old(PATH_DF::String,Nx::Int64,Ny::Int64)
    a_l = 3.335 # 2D lattice square parameter a_l = b_l = 3.335 Å = 3.335*10^(-10) m
    e_fermi = 10.3786 # Fermi energy
    kx = range(-pi, pi, length = Nx) # Range of ky values
    ky = range(-pi, pi, length = Ny) # Range of kz values
    tmp = []
    result = []
    for kk in kx
	    tmp_row = []
	    result_row = []
	    for k in ky
            println(kk, k)
		    # Compute the Hamiltonian for the given kx, ky pair
		    H = tb_Ham_old(PATH_DF, kk, k)

		    # Perform eigenvalue decomposition
		    F = eigen(H)
		    
		    # Store the eigenvalues
		    push!(tmp_row, real(F.values))
		    # Calculate the sum of squared magnitudes for the first 36 elements of each eigenvector
                    result_row_vals = [sum(abs.(F.vectors[1:36, i]).^2) for i in 1:length(F.values)]
                    push!(result_row, result_row_vals)

	    end
	    push!(tmp, tmp_row)
            push!(result, result_row)
    end
    # Flatten the nested arrays and subtract the Fermi energy
    tmp = reduce(hcat, reduce(hcat, tmp)) .- e_fermi
    result = reduce(hcat, reduce(hcat, result))
    # Return the projected plane
    return tmp, result
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

############################################################################################################

function f(w, H00, H01, H10)
    identity_matrix = Diagonal(ones(size(H00, 1)))  # Create a 36x36 identity matrix
    A = w * identity_matrix .- H00  # (w - H00)
    
    # Instead of explicitly inverting A, use \ for solving linear systems
    return H10 * (A \ H01)  # This solves A * X = H01 efficiently
end
# # Define Hamiltonian matrices (dummy matrices for example)
# H00 = rand(36, 36)  # Replace with your actual H00 matrix
# H01 = rand(36, 36)  # Replace with your actual H01 matrix
# H10 = H01'  # H10 is typically the conjugate transpose of H01

function T_int_bulk(PATH_DF::String, kx::Float64, ky::Float64)
    Matrix2DCache = Dict{String, Matrix{ComplexF64}}()
    num_2d = 36
    load_dataframes(PATH_DF)
    num_wan = maximum(
        DataFramesCache[collect(keys(DataFramesCache))[1]].n₀
        ) # Number of Wannier functions

    num_planes = num_wan ÷ num_2d
    k = 5 # 9
    Tₖ      = zeros(Complex{Float64}, num_2d, num_2d)
    d_orb0  = (k-1)*36+1:(k-1)*36+20;        #d orbitals
    p_orb0  = (k-1)*36+1+20:(k-1)*36+20+16;  #p orbitals
    m_elem0 = vcat([d_orb0, p_orb0]...)
    d_orb1  = k*36+1:k*36+20;        #d orbitals
    p_orb1  = k*36+1+20:k*36+20+16;  #p orbitals
    m_elem1 = vcat([d_orb1, p_orb1]...)
    for (filename, df) in DataFramesCache
        for i in 1:num_2d, j in 1:num_2d
            tmp = filter([:n₀, :nᵢ, :iEn] => (x, y, z) -> x == m_elem0[i] && y == m_elem1[j] && z != 0, df)
            Tₖ[i, j] += exp(im*((df.n₀[1]+0.25)*kx + df.nᵢ[1]*ky)) * (tmp.En[1] + tmp.iEn[1]*im)
        end
    end
    Matrix2DCache["T_v"] = Tₖ

    k = 6 # 10
    Tₖ      = zeros(Complex{Float64}, num_2d, num_2d)
    d_orb0  = (k-1)*36+1:(k-1)*36+20;        #d orbitals
    p_orb0  = (k-1)*36+1+20:(k-1)*36+20+16;  #p orbitals
    m_elem0 = vcat([d_orb0, p_orb0]...)
    d_orb1  = k*36+1:k*36+20;        #d orbitals
    p_orb1  = k*36+1+20:k*36+20+16;  #p orbitals
    m_elem1 = vcat([d_orb1, p_orb1]...)
    for (filename, df) in DataFramesCache
        for i in 1:num_2d, j in 1:num_2d
            tmp = filter([:n₀, :nᵢ, :iEn] => (x, y, z) -> x == m_elem0[i] && y == m_elem1[j] && z != 0, df)
            Tₖ[i, j] += exp(im*(df.n₀[1]*kx + (df.nᵢ[1]+0.25)*ky)) * (tmp.En[1] + tmp.iEn[1]*im)
        end
    end
    Matrix2DCache["T_w"] = Tₖ

    return Matrix2DCache # Returns the neares neighbor interactions from the d orbitals to the p orbitals
end
