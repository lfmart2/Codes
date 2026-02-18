using Plots
using DelimitedFiles
using HDF5
using Plots.PlotMeasures
using DataFrames
using LinearAlgebra
using LaTeXStrings
using StaticArrays

const FILE_DATA_nSOC      = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/Collinear/00_Data"
const FILE_DATA_nSOC_H    = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/Collinear_H/00_Data"
const FILE_DATA_SOC       = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/nCollinear/00_Data"
const FILE_DATA_SOC_H     = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/nCollinear_H/00_Data"
const FILE_DATA_SOC_H_S   = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell_217/nCollinear_H/00_Data"


"""
    HRData

Container for the contents of a `seedname_hr.dat` file.

Fields
-------
  num_wann   :: Int                      # number of Wannier functions
  nrpts      :: Int                      # number of distinct R–vectors
  R          :: Vector{SVector{3,Int}}   # lattice vectors (cartesian integers)
  weight     :: Vector{Int}              # degeneracy of each R
  H_R        :: Array{ComplexF64,3}      # size (num_wann,num_wann,nrpts)
"""
struct HRData
    num_wann :: Int
    nrpts    :: Int
    R        :: Matrix{Int}        # size (3, nrpts)
    weight   :: Vector{Int}
    H_R      :: Array{ComplexF64,3}  # (nw, nw, nrpts)
end

"Read seedname_hr.dat → HRData (fast, no linear scans)."
function load_hr(path::String)
    open(path, "r") do io
        readline(io)                        # header
        nw  = parse(Int, readline(io))
        nr  = parse(Int, readline(io))
        w   = Int[]
        while length(w) < nr
            append!(w, parse.(Int, split(strip(readline(io)))))
        end
        R     = zeros(Int, 3, nr)
        H_R   = Array{ComplexF64}(undef, nw, nw, nr)
        mapR  = Dict{NTuple{3,Int},Int}()
        nexti = 1
        # nr * nw * nw lines follow
        for _ in 1:(nr * nw * nw)
            vals = split(readline(io))
            Rx,Ry,Rz = parse.(Int, vals[1:3])
            ii,jj    = parse.(Int, vals[4:5])
            Re,Im    = parse.(Float64, vals[6]), parse(Float64, vals[7])
            key = (Rx,Ry,Rz)
            idx = get!(mapR, key) do
                R[:,nexti] .= (Rx,Ry,Rz)
                nexti += 1
                nexti - 1
            end
            @inbounds H_R[ii, jj, idx] = complex(Re, Im)
        end
        return HRData(nw, nr, R, w, H_R)
    end
end

"Bloch Hamiltonian H(k) = Σ_R H(R) e^{i2πk·R} / weight(R)."
function bloch_hamiltonian(hr::HRData, k::NTuple{3,Float64})
    nw, nr = hr.num_wann, hr.nrpts
    Hk = zeros(ComplexF64, nw, nw)
    @inbounds @simd for idx in 1:nr
        phase = exp(im * 2π * (k[1]*hr.R[1,idx] + k[2]*hr.R[2,idx] + k[3]*hr.R[3,idx])) / hr.weight[idx]
        Hk .+= hr.H_R[:,:,idx] * phase
    end
    return Hermitian(Hk)  # enforce numerical Hermiticity
end

#############################################################################
####### Processing the TB data for Collinear Spin without Hydrogen ##########
#############################################################################
function tb_proj_bands_nSOC(output_filename::String = joinpath(FILE_DATA_nSOC,"BANDS.h5") )
    hr_up = load_hr(joinpath(FILE_DATA_nSOC, "wannier90.up_hr.dat"))
    hr_dn = load_hr(joinpath(FILE_DATA_nSOC, "wannier90.dn_hr.dat"))
    kx = vcat(zeros(100), range(0, 0.5, length = 100), zeros(100) .+ 0.5, range(0.5, 0, length = 140))  # Range of kx values
    ky = vcat(range(0.5, 0, length = 100), zeros(100),range(0, 0.5, length = 100),range(0.5, 0, length = 140)) # Range of ky values

    k = 252
    N_tot = length(kx)

    # Compute the Hamiltonian for the given kx, ky pair

    eigvals_up, eigenvecs_up = eigen(bloch_hamiltonian(hr_up, SVector(kx[1], ky[1], 0.0)));
    eigvals_dn, eigenvecs_dn = eigen(bloch_hamiltonian(hr_dn, SVector(kx[1], ky[1], 0.0)));

    # eigvals, sorted_indices = closest_elements(real(eigvals),  e_fermi, k)
    # eigenvecs = eigenvecs[:, sorted_indices]

    # Store the eigenvalues    
    Energies_up = sort(real.(eigvals_up))
    Energies_dn = sort(real.(eigvals_dn))

    proj1_dz2_up = [abs.(eigenvecs_up[1, i]).^2 for i in 1:k]
    proj1_dz2_dn = [abs.(eigenvecs_dn[1, i]).^2 for i in 1:k]
    proj1_dxz_up = [abs.(eigenvecs_up[2, i]).^2 for i in 1:k]
    proj1_dxz_dn = [abs.(eigenvecs_dn[2, i]).^2 for i in 1:k]
    proj1_dyz_up = [abs.(eigenvecs_up[3, i]).^2 for i in 1:k]
    proj1_dyz_dn = [abs.(eigenvecs_dn[3, i]).^2 for i in 1:k]
    proj1_x2y2_up = [abs.(eigenvecs_up[4, i]).^2 for i in 1:k]
    proj1_x2y2_dn = [abs.(eigenvecs_dn[4, i]).^2 for i in 1:k]
    proj1_dxy_up = [abs.(eigenvecs_up[5, i]).^2 for i in 1:k]
    proj1_dxy_dn = [abs.(eigenvecs_dn[5, i]).^2 for i in 1:k]
    proj2_sp1_up = [abs.(eigenvecs_up[141, i]).^2 for i in 1:k]
    proj2_sp1_dn = [abs.(eigenvecs_dn[141, i]).^2 for i in 1:k]
    proj2_sp2_up = [abs.(eigenvecs_up[142, i]).^2 for i in 1:k]
    proj2_sp2_dn = [abs.(eigenvecs_dn[142, i]).^2 for i in 1:k]
    proj2_sp3_up = [abs.(eigenvecs_up[143, i]).^2 for i in 1:k]
    proj2_sp3_dn = [abs.(eigenvecs_dn[143, i]).^2 for i in 1:k]
    proj2_sp4_up = [abs.(eigenvecs_up[144, i]).^2 for i in 1:k]
    proj2_sp4_dn = [abs.(eigenvecs_dn[144, i]).^2 for i in 1:k]


    Projections1_up = [sum(abs.(eigenvecs_up[1:5, i]).^2) for i in 1:k]
    Projections2_up = [sum(abs.(eigenvecs_up[141:144, i]).^2) for i in 1:k]
    Projections3_up = [sum(abs.(eigenvecs_up[6:10, i]).^2) for i in 1:k]
    Projections4_up = [sum(abs.(eigenvecs_up[145:148, i]).^2) for i in 1:k]
    Projections5_up = [sum(abs.(eigenvecs_up[11:15, i]).^2) for i in 1:k]
    Projections6_up = [sum(abs.(eigenvecs_up[149:152, i]).^2) for i in 1:k]
    Projections7_up = [sum(abs.(eigenvecs_up[16:20, i]).^2) for i in 1:k]
    Projections8_up = [sum(abs.(eigenvecs_up[153:156, i]).^2) for i in 1:k]

    Projections1_dn = [sum(abs.(eigenvecs_dn[1:5, i]).^2) for i in 1:k]
    Projections2_dn = [sum(abs.(eigenvecs_dn[141:144, i]).^2) for i in 1:k]
    Projections3_dn = [sum(abs.(eigenvecs_dn[6:10, i]).^2) for i in 1:k]
    Projections4_dn = [sum(abs.(eigenvecs_dn[145:148, i]).^2) for i in 1:k]
    Projections5_dn = [sum(abs.(eigenvecs_dn[11:15, i]).^2) for i in 1:k]
    Projections6_dn = [sum(abs.(eigenvecs_dn[149:152, i]).^2) for i in 1:k]
    Projections7_dn = [sum(abs.(eigenvecs_dn[16:20, i]).^2) for i in 1:k]
    Projections8_dn = [sum(abs.(eigenvecs_dn[153:156, i]).^2) for i in 1:k]
    for i in 2:N_tot
        # Perform eigenvalue decomposition
        eigvals_up, eigenvecs_up = eigen(bloch_hamiltonian(hr_up, SVector(kx[i], ky[i], 0.0)));
        eigvals_dn, eigenvecs_dn = eigen(bloch_hamiltonian(hr_dn, SVector(kx[i], ky[i], 0.0)));        # eigvals, sorted_indices = closest_elements(real(eigvals),  e_fermi, k)

        # Store the eigenvalues
        Energies_up = hcat(Energies_up, sort(real.(eigvals_up)))
        Energies_dn = hcat(Energies_dn, sort(real.(eigvals_dn)))

        proj1_dz2_up = hcat(proj1_dz2_up, [abs.(eigenvecs_up[1, i]).^2 for i in 1:k])
        proj1_dz2_dn = hcat(proj1_dz2_dn, [abs.(eigenvecs_dn[1, i]).^2 for i in 1:k])
        proj1_dxz_up = hcat(proj1_dxz_up, [abs.(eigenvecs_up[2, i]).^2 for i in 1:k])
        proj1_dxz_dn = hcat(proj1_dxz_dn, [abs.(eigenvecs_dn[2, i]).^2 for i in 1:k])
        proj1_dyz_up = hcat(proj1_dyz_up, [abs.(eigenvecs_up[3, i]).^2 for i in 1:k])
        proj1_dyz_dn = hcat(proj1_dyz_dn, [abs.(eigenvecs_dn[3, i]).^2 for i in 1:k])
        proj1_x2y2_up = hcat(proj1_x2y2_up, [abs.(eigenvecs_up[4, i]).^2 for i in 1:k])
        proj1_x2y2_dn = hcat(proj1_x2y2_dn, [abs.(eigenvecs_dn[4, i]).^2 for i in 1:k])
        proj1_dxy_up = hcat(proj1_dxy_up, [abs.(eigenvecs_up[5, i]).^2 for i in 1:k])
        proj1_dxy_dn = hcat(proj1_dxy_dn, [abs.(eigenvecs_dn[5, i]).^2 for i in 1:k])
        proj2_sp1_up = hcat(proj2_sp1_up, [abs.(eigenvecs_up[141, i]).^2 for i in 1:k])
        proj2_sp1_dn = hcat(proj2_sp1_dn, [abs.(eigenvecs_dn[141, i]).^2 for i in 1:k])
        proj2_sp2_up = hcat(proj2_sp2_up, [abs.(eigenvecs_up[142, i]).^2 for i in 1:k])
        proj2_sp2_dn = hcat(proj2_sp2_dn, [abs.(eigenvecs_dn[142, i]).^2 for i in 1:k])
        proj2_sp3_up = hcat(proj2_sp3_up, [abs.(eigenvecs_up[143, i]).^2 for i in 1:k])
        proj2_sp3_dn = hcat(proj2_sp3_dn, [abs.(eigenvecs_dn[143, i]).^2 for i in 1:k])
        proj2_sp4_up = hcat(proj2_sp4_up, [abs.(eigenvecs_up[144, i]).^2 for i in 1:k])
        proj2_sp4_dn = hcat(proj2_sp4_dn, [abs.(eigenvecs_dn[144, i]).^2 for i in 1:k])

        Projections1_up = hcat(Projections1_up, [sum(abs.(eigenvecs_up[1:5, i]).^2) for i in 1:k])
        Projections2_up = hcat(Projections2_up, [sum(abs.(eigenvecs_up[141:144, i]).^2) for i in 1:k])
        Projections3_up = hcat(Projections3_up, [sum(abs.(eigenvecs_up[6:10, i]).^2) for i in 1:k])
        Projections4_up = hcat(Projections4_up, [sum(abs.(eigenvecs_up[145:148, i]).^2) for i in 1:k])
        Projections5_up = hcat(Projections5_up, [sum(abs.(eigenvecs_up[11:15, i]).^2) for i in 1:k])
        Projections6_up = hcat(Projections6_up, [sum(abs.(eigenvecs_up[149:152, i]).^2) for i in 1:k])
        Projections7_up = hcat(Projections7_up, [sum(abs.(eigenvecs_up[16:20, i]).^2) for i in 1:k])
        Projections8_up = hcat(Projections8_up, [sum(abs.(eigenvecs_up[153:156, i]).^2) for i in 1:k])

        Projections1_dn = hcat(Projections1_dn, [sum(abs.(eigenvecs_dn[1:5, i]).^2) for i in 1:k])
        Projections2_dn = hcat(Projections2_dn, [sum(abs.(eigenvecs_dn[141:144, i]).^2) for i in 1:k])
        Projections3_dn = hcat(Projections3_dn, [sum(abs.(eigenvecs_dn[6:10, i]).^2) for i in 1:k])
        Projections4_dn = hcat(Projections4_dn, [sum(abs.(eigenvecs_dn[145:148, i]).^2) for i in 1:k])
        Projections5_dn = hcat(Projections5_dn, [sum(abs.(eigenvecs_dn[11:15, i]).^2) for i in 1:k])
        Projections6_dn = hcat(Projections6_dn, [sum(abs.(eigenvecs_dn[149:152, i]).^2) for i in 1:k])
        Projections7_dn = hcat(Projections7_dn, [sum(abs.(eigenvecs_dn[16:20, i]).^2) for i in 1:k])
        Projections8_dn = hcat(Projections8_dn, [sum(abs.(eigenvecs_dn[153:156, i]).^2) for i in 1:k])
    end

    h5open(output_filename, "w") do file
        file["eigenvalues_up"] = Energies_up
        file["eigenvalues_dn"] = Energies_dn

        file["proj1_dxy_up"] = proj1_dxy_up
        file["proj1_dxy_dn"] = proj1_dxy_dn
        file["proj1_dyz_up"] = proj1_dyz_up
        file["proj1_dyz_dn"] = proj1_dyz_dn
        file["proj1_dz2_up"] = proj1_dz2_up
        file["proj1_dz2_dn"] = proj1_dz2_dn
        file["proj1_dxz_up"] = proj1_dxz_up
        file["proj1_dxz_dn"] = proj1_dxz_dn
        file["proj1_x2y2_up"] = proj1_x2y2_up
        file["proj1_x2y2_dn"] = proj1_x2y2_dn
        file["proj2_sp1_up"] = proj2_sp1_up
        file["proj2_sp1_dn"] = proj2_sp1_dn
        file["proj2_sp2_up"] = proj2_sp2_up
        file["proj2_sp2_dn"] = proj2_sp2_dn
        file["proj2_sp3_up"] = proj2_sp3_up
        file["proj2_sp3_dn"] = proj2_sp3_dn
        file["proj2_sp4_up"] = proj2_sp4_up
        file["proj2_sp4_dn"] = proj2_sp4_dn

        file["proj1_up"] = Projections1_up
        file["proj2_up"] = Projections2_up
        file["proj3_up"] = Projections3_up
        file["proj4_up"] = Projections4_up
        file["proj5_up"] = Projections5_up
        file["proj6_up"] = Projections6_up
        file["proj7_up"] = Projections7_up
        file["proj8_up"] = Projections8_up

        file["proj1_dn"] = Projections1_dn
        file["proj2_dn"] = Projections2_dn
        file["proj3_dn"] = Projections3_dn
        file["proj4_dn"] = Projections4_dn
        file["proj5_dn"] = Projections5_dn
        file["proj6_dn"] = Projections6_dn
        file["proj7_dn"] = Projections7_dn
        file["proj8_dn"] = Projections8_dn
    end

end

#############################################################################
######### Processing the TB data for Collinear Spin with Hydrogen ###########
#############################################################################
function tb_proj_bands_nSOC_H(output_filename::String = joinpath(FILE_DATA_nSOC_H,"BANDS_H.h5") )
    hr_H_up = load_hr(joinpath(FILE_DATA_nSOC_H, "wannier90.up_hr.dat"))
    hr_H_dn = load_hr(joinpath(FILE_DATA_nSOC_H, "wannier90.dn_hr.dat"))
    kx = vcat(zeros(100), range(0, 0.5, length = 100), zeros(100) .+ 0.5, range(0.5, 0, length = 140))  # Range of kx values
    ky = vcat(range(0.5, 0, length = 100), zeros(100),range(0, 0.5, length = 100),range(0.5, 0, length = 140)) # Range of ky values

    k = 253
    N_tot = length(kx)

    # Compute the Hamiltonian for the given kx, ky pair

    eigvals_up, eigenvecs_up = eigen(bloch_hamiltonian(hr_H_up, SVector(kx[1], ky[1], 0.0)));
    eigvals_dn, eigenvecs_dn = eigen(bloch_hamiltonian(hr_H_dn, SVector(kx[1], ky[1], 0.0)));

    # eigvals, sorted_indices = closest_elements(real(eigvals),  e_fermi, k)
    # eigenvecs = eigenvecs[:, sorted_indices]

    # Store the eigenvalues    
    Energies_up = sort(real.(eigvals_up))
    Energies_dn = sort(real.(eigvals_dn))

    proj0_H_up = [abs.(eigenvecs_up[1,i]).^2 for i in 1:k]
	proj0_H_dn = [abs.(eigenvecs_dn[1,i]).^2 for i in 1:k]
    proj1_dz2_up = [abs.(eigenvecs_up[2, i]).^2 for i in 1:k]
    proj1_dz2_dn = [abs.(eigenvecs_dn[2, i]).^2 for i in 1:k]
    proj1_dxz_up = [abs.(eigenvecs_up[3, i]).^2 for i in 1:k]
    proj1_dxz_dn = [abs.(eigenvecs_dn[3, i]).^2 for i in 1:k]
    proj1_dyz_up = [abs.(eigenvecs_up[4, i]).^2 for i in 1:k]
    proj1_dyz_dn = [abs.(eigenvecs_dn[4, i]).^2 for i in 1:k]
    proj1_x2y2_up = [abs.(eigenvecs_up[5, i]).^2 for i in 1:k]
    proj1_x2y2_dn = [abs.(eigenvecs_dn[5, i]).^2 for i in 1:k]
    proj1_dxy_up = [abs.(eigenvecs_up[6, i]).^2 for i in 1:k]
    proj1_dxy_dn = [abs.(eigenvecs_dn[6, i]).^2 for i in 1:k]
    proj2_sp1_up = [abs.(eigenvecs_up[142, i]).^2 for i in 1:k]
    proj2_sp1_dn = [abs.(eigenvecs_dn[142, i]).^2 for i in 1:k]
    proj2_sp2_up = [abs.(eigenvecs_up[143, i]).^2 for i in 1:k]
    proj2_sp2_dn = [abs.(eigenvecs_dn[143, i]).^2 for i in 1:k]
    proj2_sp3_up = [abs.(eigenvecs_up[144, i]).^2 for i in 1:k]
    proj2_sp3_dn = [abs.(eigenvecs_dn[144, i]).^2 for i in 1:k]
    proj2_sp4_up = [abs.(eigenvecs_up[145, i]).^2 for i in 1:k]
    proj2_sp4_dn = [abs.(eigenvecs_dn[145, i]).^2 for i in 1:k]


    Projections1_up = [sum(abs.(eigenvecs_up[2:6, i]).^2) for i in 1:k]
    Projections2_up = [sum(abs.(eigenvecs_up[142:145, i]).^2) for i in 1:k]
    Projections3_up = [sum(abs.(eigenvecs_up[7:11, i]).^2) for i in 1:k]
    Projections4_up = [sum(abs.(eigenvecs_up[146:149, i]).^2) for i in 1:k]
    Projections5_up = [sum(abs.(eigenvecs_up[12:16, i]).^2) for i in 1:k]
    Projections6_up = [sum(abs.(eigenvecs_up[150:153, i]).^2) for i in 1:k]
    Projections7_up = [sum(abs.(eigenvecs_up[17:21, i]).^2) for i in 1:k]
    Projections8_up = [sum(abs.(eigenvecs_up[154:157, i]).^2) for i in 1:k]

    Projections1_dn = [sum(abs.(eigenvecs_dn[2:6, i]).^2) for i in 1:k]
    Projections2_dn = [sum(abs.(eigenvecs_dn[142:145, i]).^2) for i in 1:k]
    Projections3_dn = [sum(abs.(eigenvecs_dn[7:11, i]).^2) for i in 1:k]
    Projections4_dn = [sum(abs.(eigenvecs_dn[146:149, i]).^2) for i in 1:k]
    Projections5_dn = [sum(abs.(eigenvecs_dn[12:16, i]).^2) for i in 1:k]
    Projections6_dn = [sum(abs.(eigenvecs_dn[150:153, i]).^2) for i in 1:k]
    Projections7_dn = [sum(abs.(eigenvecs_dn[17:21, i]).^2) for i in 1:k]
    Projections8_dn = [sum(abs.(eigenvecs_dn[154:157, i]).^2) for i in 1:k]
    for i in 2:length(kx)
        # Perform eigenvalue decomposition
        eigvals_up, eigenvecs_up = eigen(bloch_hamiltonian(hr_H_up, SVector(kx[i], ky[i], 0.0)));
        eigvals_dn, eigenvecs_dn = eigen(bloch_hamiltonian(hr_H_dn, SVector(kx[i], ky[i], 0.0)));        # eigvals, sorted_indices = closest_elements(real(eigvals),  e_fermi, k)

        # Store the eigenvalues
        Energies_up = hcat(Energies_up, sort(real.(eigvals_up)))
        Energies_dn = hcat(Energies_dn, sort(real.(eigvals_dn)))

        proj0_H_up = hcat(proj0_H_up, [abs.(eigenvecs_up[1,i]).^2 for i in 1:k])
        proj0_H_dn = hcat(proj0_H_dn, [abs.(eigenvecs_dn[1,i]).^2 for i in 1:k])
        proj1_dz2_up = hcat(proj1_dz2_up, [abs.(eigenvecs_up[2, i]).^2 for i in 1:k])
        proj1_dz2_dn = hcat(proj1_dz2_dn, [abs.(eigenvecs_dn[2, i]).^2 for i in 1:k])
        proj1_dxz_up = hcat(proj1_dxz_up, [abs.(eigenvecs_up[3, i]).^2 for i in 1:k])
        proj1_dxz_dn = hcat(proj1_dxz_dn, [abs.(eigenvecs_dn[3, i]).^2 for i in 1:k])
        proj1_dyz_up = hcat(proj1_dyz_up, [abs.(eigenvecs_up[4, i]).^2 for i in 1:k])
        proj1_dyz_dn = hcat(proj1_dyz_dn, [abs.(eigenvecs_dn[4, i]).^2 for i in 1:k])
        proj1_x2y2_up = hcat(proj1_x2y2_up, [abs.(eigenvecs_up[5, i]).^2 for i in 1:k])
        proj1_x2y2_dn = hcat(proj1_x2y2_dn, [abs.(eigenvecs_dn[5, i]).^2 for i in 1:k])
        proj1_dxy_up = hcat(proj1_dxy_up, [abs.(eigenvecs_up[6, i]).^2 for i in 1:k])
        proj1_dxy_dn = hcat(proj1_dxy_dn, [abs.(eigenvecs_dn[7, i]).^2 for i in 1:k])
        proj2_sp1_up = hcat(proj2_sp1_up, [abs.(eigenvecs_up[142, i]).^2 for i in 1:k])
        proj2_sp1_dn = hcat(proj2_sp1_dn, [abs.(eigenvecs_dn[142, i]).^2 for i in 1:k])
        proj2_sp2_up = hcat(proj2_sp2_up, [abs.(eigenvecs_up[143, i]).^2 for i in 1:k])
        proj2_sp2_dn = hcat(proj2_sp2_dn, [abs.(eigenvecs_dn[143, i]).^2 for i in 1:k])
        proj2_sp3_up = hcat(proj2_sp3_up, [abs.(eigenvecs_up[144, i]).^2 for i in 1:k])
        proj2_sp3_dn = hcat(proj2_sp3_dn, [abs.(eigenvecs_dn[144, i]).^2 for i in 1:k])
        proj2_sp4_up = hcat(proj2_sp4_up, [abs.(eigenvecs_up[145, i]).^2 for i in 1:k])
        proj2_sp4_dn = hcat(proj2_sp4_dn, [abs.(eigenvecs_dn[145, i]).^2 for i in 1:k])

        Projections1_up = hcat(Projections1_up, [sum(abs.(eigenvecs_up[2:6, i]).^2) for i in 1:k])
        Projections2_up = hcat(Projections2_up, [sum(abs.(eigenvecs_up[142:145, i]).^2) for i in 1:k])
        Projections3_up = hcat(Projections3_up, [sum(abs.(eigenvecs_up[7:11, i]).^2) for i in 1:k])
        Projections4_up = hcat(Projections4_up, [sum(abs.(eigenvecs_up[146:149, i]).^2) for i in 1:k])
        Projections5_up = hcat(Projections5_up, [sum(abs.(eigenvecs_up[12:16, i]).^2) for i in 1:k])
        Projections6_up = hcat(Projections6_up, [sum(abs.(eigenvecs_up[150:153, i]).^2) for i in 1:k])
        Projections7_up = hcat(Projections7_up, [sum(abs.(eigenvecs_up[17:21, i]).^2) for i in 1:k])
        Projections8_up = hcat(Projections8_up, [sum(abs.(eigenvecs_up[154:157, i]).^2) for i in 1:k])

        Projections1_dn = hcat(Projections1_dn, [sum(abs.(eigenvecs_dn[2:6, i]).^2) for i in 1:k])
        Projections2_dn = hcat(Projections2_dn, [sum(abs.(eigenvecs_dn[142:145, i]).^2) for i in 1:k])
        Projections3_dn = hcat(Projections3_dn, [sum(abs.(eigenvecs_dn[7:11, i]).^2) for i in 1:k])
        Projections4_dn = hcat(Projections4_dn, [sum(abs.(eigenvecs_dn[146:149, i]).^2) for i in 1:k])
        Projections5_dn = hcat(Projections5_dn, [sum(abs.(eigenvecs_dn[12:16, i]).^2) for i in 1:k])
        Projections6_dn = hcat(Projections6_dn, [sum(abs.(eigenvecs_dn[150:153, i]).^2) for i in 1:k])
        Projections7_dn = hcat(Projections7_dn, [sum(abs.(eigenvecs_dn[17:21, i]).^2) for i in 1:k])
        Projections8_dn = hcat(Projections8_dn, [sum(abs.(eigenvecs_dn[154:157, i]).^2) for i in 1:k])
    end

    h5open(output_filename, "w") do file
        file["eigenvalues_up"] = Energies_up
        file["eigenvalues_dn"] = Energies_dn

        file["proj_H_up"] = proj0_H_up
		file["proj_H_dn"] = proj0_H_dn
        file["proj1_dxy_up"] = proj1_dxy_up
        file["proj1_dxy_dn"] = proj1_dxy_dn
        file["proj1_dyz_up"] = proj1_dyz_up
        file["proj1_dyz_dn"] = proj1_dyz_dn
        file["proj1_dz2_up"] = proj1_dz2_up
        file["proj1_dz2_dn"] = proj1_dz2_dn
        file["proj1_dxz_up"] = proj1_dxz_up
        file["proj1_dxz_dn"] = proj1_dxz_dn
        file["proj1_x2y2_up"] = proj1_x2y2_up
        file["proj1_x2y2_dn"] = proj1_x2y2_dn
        file["proj2_sp1_up"] = proj2_sp1_up
        file["proj2_sp1_dn"] = proj2_sp1_dn
        file["proj2_sp2_up"] = proj2_sp2_up
        file["proj2_sp2_dn"] = proj2_sp2_dn
        file["proj2_sp3_up"] = proj2_sp3_up
        file["proj2_sp3_dn"] = proj2_sp3_dn
        file["proj2_sp4_up"] = proj2_sp4_up
        file["proj2_sp4_dn"] = proj2_sp4_dn

        file["proj1_up"] = Projections1_up
        file["proj2_up"] = Projections2_up
        file["proj3_up"] = Projections3_up
        file["proj4_up"] = Projections4_up
        file["proj5_up"] = Projections5_up
        file["proj6_up"] = Projections6_up
        file["proj7_up"] = Projections7_up
        file["proj8_up"] = Projections8_up

        file["proj1_dn"] = Projections1_dn
        file["proj2_dn"] = Projections2_dn
        file["proj3_dn"] = Projections3_dn
        file["proj4_dn"] = Projections4_dn
        file["proj5_dn"] = Projections5_dn
        file["proj6_dn"] = Projections6_dn
        file["proj7_dn"] = Projections7_dn
        file["proj8_dn"] = Projections8_dn
    end

end

#############################################################################
######## Processing the TB data for nonCollinear Spin wihtout Hydrogen #########
#############################################################################
function tb_proj_bands_SOC(output_filename::String = joinpath(@__DIR__,"BANDS_gpu.h5") )
    hr = load_hr(joinpath(FILE_DATA_SOC,"wannier90_hr.dat")  )
    # kx   = vcat(zeros(100), range(0, 0.5, length = 100), zeros(100) .+ 0.5, range(0.5, 0, length = 140))  # Range of kx values
    # ky   = vcat(range(0.5, 0, length = 100), zeros(100),range(0, 0.5, length = 100),range(0.5, 0, length = 140)) # Range of ky values

    kx, ky = tst()
    k = 504
    N_tot = length(kx)

    # Compute the Hamiltonian for the given kx, ky pair

    eigvals, eigenvecs = eigen(bloch_hamiltonian(hr, (kx[1], ky[1], 0.0)));

    # Store the eigenvalues    
    Energies = sort(real.(eigvals))

    proj1_dz2_up = [abs.(eigenvecs[1, i]).^2 for i in 1:k]
    proj1_dz2_dn = [abs.(eigenvecs[2, i]).^2 for i in 1:k]
    proj1_dxz_up = [abs.(eigenvecs[3, i]).^2 for i in 1:k]
    proj1_dxz_dn = [abs.(eigenvecs[4, i]).^2 for i in 1:k]
    proj1_dyz_up = [abs.(eigenvecs[5, i]).^2 for i in 1:k]
    proj1_dyz_dn = [abs.(eigenvecs[6, i]).^2 for i in 1:k]
    proj1_x2y2_up = [abs.(eigenvecs[7, i]).^2 for i in 1:k]
    proj1_x2y2_dn = [abs.(eigenvecs[8, i]).^2 for i in 1:k]
    proj1_dxy_up = [abs.(eigenvecs[9, i]).^2 for i in 1:k]
    proj1_dxy_dn = [abs.(eigenvecs[10, i]).^2 for i in 1:k]
    proj2_sp1_up = [abs.(eigenvecs[281, i]).^2 for i in 1:k]
    proj2_sp1_dn = [abs.(eigenvecs[282, i]).^2 for i in 1:k]
    proj2_sp2_up = [abs.(eigenvecs[283, i]).^2 for i in 1:k]
    proj2_sp2_dn = [abs.(eigenvecs[284, i]).^2 for i in 1:k]
    proj2_sp3_up = [abs.(eigenvecs[285, i]).^2 for i in 1:k]
    proj2_sp3_dn = [abs.(eigenvecs[286, i]).^2 for i in 1:k]
    proj2_sp4_up = [abs.(eigenvecs[287, i]).^2 for i in 1:k]
    proj2_sp4_dn = [abs.(eigenvecs[288, i]).^2 for i in 1:k]


    Projections1 = [sum(abs.(eigenvecs[1:10, i]).^2) for i in 1:k]
    Projections2 = [sum(abs.(eigenvecs[281:288, i]).^2) for i in 1:k]
    Projections3 = [sum(abs.(eigenvecs[11:20, i]).^2) for i in 1:k]
    Projections4 = [sum(abs.(eigenvecs[289:296, i]).^2) for i in 1:k]
    Projections5 = [sum(abs.(eigenvecs[21:30, i]).^2) for i in 1:k]
    Projections6 = [sum(abs.(eigenvecs[297:304, i]).^2) for i in 1:k]
    Projections7 = [sum(abs.(eigenvecs[31:40, i]).^2) for i in 1:k]
    Projections8 = [sum(abs.(eigenvecs[305:312, i]).^2) for i in 1:k]
    for i in 2:N_tot
        # Perform eigenvalue decomposition
        eigvals, eigenvecs = eigen(bloch_hamiltonian(hr, (kx[i], ky[i], 0.0)));

        # Store the eigenvalues
        Energies = hcat(Energies, sort(real.(eigvals)))

        proj1_dz2_up = hcat(proj1_dz2_up, [abs.(eigenvecs[1, i]).^2 for i in 1:k])
        proj1_dz2_dn = hcat(proj1_dz2_dn, [abs.(eigenvecs[2, i]).^2 for i in 1:k])
        proj1_dxz_up = hcat(proj1_dxz_up, [abs.(eigenvecs[3, i]).^2 for i in 1:k])
        proj1_dxz_dn = hcat(proj1_dxz_dn, [abs.(eigenvecs[4, i]).^2 for i in 1:k])
        proj1_dyz_up = hcat(proj1_dyz_up, [abs.(eigenvecs[5, i]).^2 for i in 1:k])
        proj1_dyz_dn = hcat(proj1_dyz_dn, [abs.(eigenvecs[6, i]).^2 for i in 1:k])
        proj1_x2y2_up = hcat(proj1_x2y2_up, [abs.(eigenvecs[7, i]).^2 for i in 1:k])
        proj1_x2y2_dn = hcat(proj1_x2y2_dn, [abs.(eigenvecs[8, i]).^2 for i in 1:k])
        proj1_dxy_up = hcat(proj1_dxy_up, [abs.(eigenvecs[9, i]).^2 for i in 1:k])
        proj1_dxy_dn = hcat(proj1_dxy_dn, [abs.(eigenvecs[10, i]).^2 for i in 1:k])
        proj2_sp1_up = hcat(proj2_sp1_up, [abs.(eigenvecs[281, i]).^2 for i in 1:k])
        proj2_sp1_dn = hcat(proj2_sp1_dn, [abs.(eigenvecs[282, i]).^2 for i in 1:k])
        proj2_sp2_up = hcat(proj2_sp2_up, [abs.(eigenvecs[283, i]).^2 for i in 1:k])
        proj2_sp2_dn = hcat(proj2_sp2_dn, [abs.(eigenvecs[284, i]).^2 for i in 1:k])
        proj2_sp3_up = hcat(proj2_sp3_up, [abs.(eigenvecs[285, i]).^2 for i in 1:k])
        proj2_sp3_dn = hcat(proj2_sp3_dn, [abs.(eigenvecs[286, i]).^2 for i in 1:k])
        proj2_sp4_up = hcat(proj2_sp4_up, [abs.(eigenvecs[287, i]).^2 for i in 1:k])
        proj2_sp4_dn = hcat(proj2_sp4_dn, [abs.(eigenvecs[288, i]).^2 for i in 1:k])

        Projections1 = hcat(Projections1, [sum(abs.(eigenvecs[1:10, i]).^2) for i in 1:k])
        Projections2 = hcat(Projections2, [sum(abs.(eigenvecs[281:288, i]).^2) for i in 1:k])
        Projections3 = hcat(Projections3, [sum(abs.(eigenvecs[11:20, i]).^2) for i in 1:k])
        Projections4 = hcat(Projections4, [sum(abs.(eigenvecs[289:296, i]).^2) for i in 1:k])
        Projections5 = hcat(Projections5, [sum(abs.(eigenvecs[21:30, i]).^2) for i in 1:k])
        Projections6 = hcat(Projections6, [sum(abs.(eigenvecs[297:304, i]).^2) for i in 1:k])
        Projections7 = hcat(Projections7, [sum(abs.(eigenvecs[33:40, i]).^2) for i in 1:k])
        Projections8 = hcat(Projections8, [sum(abs.(eigenvecs[305:312, i]).^2) for i in 1:k])
    end

    h5open(output_filename, "w") do file
        file["eigenvalues"] = Energies

        file["proj1_dxy_up"] = proj1_dxy_up
        file["proj1_dxy_dn"] = proj1_dxy_dn
        file["proj1_dyz_up"] = proj1_dyz_up
        file["proj1_dyz_dn"] = proj1_dyz_dn
        file["proj1_dz2_up"] = proj1_dz2_up
        file["proj1_dz2_dn"] = proj1_dz2_dn
        file["proj1_dxz_up"] = proj1_dxz_up
        file["proj1_dxz_dn"] = proj1_dxz_dn
        file["proj1_x2y2_up"] = proj1_x2y2_up
        file["proj1_x2y2_dn"] = proj1_x2y2_dn
        file["proj2_sp1_up"] = proj2_sp1_up
        file["proj2_sp1_dn"] = proj2_sp1_dn
        file["proj2_sp2_up"] = proj2_sp2_up
        file["proj2_sp2_dn"] = proj2_sp2_dn
        file["proj2_sp3_up"] = proj2_sp3_up
        file["proj2_sp3_dn"] = proj2_sp3_dn
        file["proj2_sp4_up"] = proj2_sp4_up
        file["proj2_sp4_dn"] = proj2_sp4_dn

        file["proj1"] = Projections1
        file["proj2"] = Projections2
        file["proj3"] = Projections3
        file["proj4"] = Projections4
        file["proj5"] = Projections5
        file["proj6"] = Projections6
        file["proj7"] = Projections7
        file["proj8"] = Projections8
    end

end

#############################################################################
########## Processing the TB data for Collinear Spin with Hydrogen ##########
#############################################################################
function tb_proj_bands_SOC_H(output_filename::String = joinpath(@__DIR__,"BANDS_H2.h5") )
    hr_H = load_hr(joinpath(FILE_DATA_SOC_H,"wannier90_hr.dat")  )
    # kx   = vcat(zeros(100), range(0, 0.5, length = 100), zeros(100) .+ 0.5, range(0.5, 0, length = 140))  # Range of kx values
    # ky   = vcat(range(0.5, 0, length = 100), zeros(100),range(0, 0.5, length = 100),range(0.5, 0, length = 140)) # Range of ky values

    kx, ky = tst()
    k = 506
    N_tot = length(kx)

    # Compute the Hamiltonian for the given kx, ky pair

    eigvals, eigenvecs = eigen(bloch_hamiltonian(hr_H, (kx[1], ky[1], 0.0)));

    # eigvals, sorted_indices = closest_elements(real(eigvals),  e_fermi, k)
    # eigenvecs = eigenvecs[:, sorted_indices]

    # Store the eigenvalues    
    Energies = sort(real.(eigvals))

    proj0_H_up = [abs.(eigenvecs[1,i]).^2 for i in 1:k]
	proj0_H_dn = [abs.(eigenvecs[2,i]).^2 for i in 1:k]
    proj1_dz2_up = [abs.(eigenvecs[3, i]).^2 for i in 1:k]
    proj1_dz2_dn = [abs.(eigenvecs[4, i]).^2 for i in 1:k]
    proj1_dxz_up = [abs.(eigenvecs[5, i]).^2 for i in 1:k]
    proj1_dxz_dn = [abs.(eigenvecs[6, i]).^2 for i in 1:k]
    proj1_dyz_up = [abs.(eigenvecs[7, i]).^2 for i in 1:k]
    proj1_dyz_dn = [abs.(eigenvecs[8, i]).^2 for i in 1:k]
    proj1_x2y2_up = [abs.(eigenvecs[9, i]).^2 for i in 1:k]
    proj1_x2y2_dn = [abs.(eigenvecs[10, i]).^2 for i in 1:k]
    proj1_dxy_up = [abs.(eigenvecs[11, i]).^2 for i in 1:k]
    proj1_dxy_dn = [abs.(eigenvecs[12, i]).^2 for i in 1:k]
    proj2_sp1_up = [abs.(eigenvecs[283, i]).^2 for i in 1:k]
    proj2_sp1_dn = [abs.(eigenvecs[284, i]).^2 for i in 1:k]
    proj2_sp2_up = [abs.(eigenvecs[285, i]).^2 for i in 1:k]
    proj2_sp2_dn = [abs.(eigenvecs[286, i]).^2 for i in 1:k]
    proj2_sp3_up = [abs.(eigenvecs[287, i]).^2 for i in 1:k]
    proj2_sp3_dn = [abs.(eigenvecs[288, i]).^2 for i in 1:k]
    proj2_sp4_up = [abs.(eigenvecs[289, i]).^2 for i in 1:k]
    proj2_sp4_dn = [abs.(eigenvecs[290, i]).^2 for i in 1:k]


    Projections1 = [sum(abs.(eigenvecs[3:12, i]).^2) for i in 1:k]
    Projections2 = [sum(abs.(eigenvecs[283:290, i]).^2) for i in 1:k]
    Projections3 = [sum(abs.(eigenvecs[13:22, i]).^2) for i in 1:k]
    Projections4 = [sum(abs.(eigenvecs[291:298, i]).^2) for i in 1:k]
    Projections5 = [sum(abs.(eigenvecs[23:32, i]).^2) for i in 1:k]
    Projections6 = [sum(abs.(eigenvecs[299:306, i]).^2) for i in 1:k]
    Projections7 = [sum(abs.(eigenvecs[33:42, i]).^2) for i in 1:k]
    Projections8 = [sum(abs.(eigenvecs[307:314, i]).^2) for i in 1:k]
    for i in 2:length(kx)
        # Perform eigenvalue decomposition
        eigvals, eigenvecs = eigen(bloch_hamiltonian(hr_H, (kx[i], ky[i], 0.0)));
        # eigvals, sorted_indices = closest_elements(real(eigvals),  e_fermi, k)
        # eigenvecs = eigenvecs[:, sorted_indices]

        # Store the eigenvalues
        Energies = hcat(Energies, sort(real.(eigvals)))

        proj0_H_up = hcat(proj0_H_up, [abs.(eigenvecs[1,i]).^2 for i in 1:k])
        proj0_H_dn = hcat(proj0_H_dn, [abs.(eigenvecs[2,i]).^2 for i in 1:k])
        proj1_dz2_up = hcat(proj1_dz2_up, [abs.(eigenvecs[3, i]).^2 for i in 1:k])
        proj1_dz2_dn = hcat(proj1_dz2_dn, [abs.(eigenvecs[4, i]).^2 for i in 1:k])
        proj1_dxz_up = hcat(proj1_dxz_up, [abs.(eigenvecs[5, i]).^2 for i in 1:k])
        proj1_dxz_dn = hcat(proj1_dxz_dn, [abs.(eigenvecs[6, i]).^2 for i in 1:k])
        proj1_dyz_up = hcat(proj1_dyz_up, [abs.(eigenvecs[7, i]).^2 for i in 1:k])
        proj1_dyz_dn = hcat(proj1_dyz_dn, [abs.(eigenvecs[8, i]).^2 for i in 1:k])
        proj1_x2y2_up = hcat(proj1_x2y2_up, [abs.(eigenvecs[9, i]).^2 for i in 1:k])
        proj1_x2y2_dn = hcat(proj1_x2y2_dn, [abs.(eigenvecs[10, i]).^2 for i in 1:k])
        proj1_dxy_up = hcat(proj1_dxy_up, [abs.(eigenvecs[11, i]).^2 for i in 1:k])
        proj1_dxy_dn = hcat(proj1_dxy_dn, [abs.(eigenvecs[12, i]).^2 for i in 1:k])
        proj2_sp1_up = hcat(proj2_sp1_up, [abs.(eigenvecs[283, i]).^2 for i in 1:k])
        proj2_sp1_dn = hcat(proj2_sp1_dn, [abs.(eigenvecs[284, i]).^2 for i in 1:k])
        proj2_sp2_up = hcat(proj2_sp2_up, [abs.(eigenvecs[285, i]).^2 for i in 1:k])
        proj2_sp2_dn = hcat(proj2_sp2_dn, [abs.(eigenvecs[286, i]).^2 for i in 1:k])
        proj2_sp3_up = hcat(proj2_sp3_up, [abs.(eigenvecs[287, i]).^2 for i in 1:k])
        proj2_sp3_dn = hcat(proj2_sp3_dn, [abs.(eigenvecs[288, i]).^2 for i in 1:k])
        proj2_sp4_up = hcat(proj2_sp4_up, [abs.(eigenvecs[289, i]).^2 for i in 1:k])
        proj2_sp4_dn = hcat(proj2_sp4_dn, [abs.(eigenvecs[290, i]).^2 for i in 1:k])

        Projections1 = hcat(Projections1, [sum(abs.(eigenvecs[3:12, i]).^2) for i in 1:k])
        Projections2 = hcat(Projections2, [sum(abs.(eigenvecs[283:290, i]).^2) for i in 1:k])
        Projections3 = hcat(Projections3, [sum(abs.(eigenvecs[13:22, i]).^2) for i in 1:k])
        Projections4 = hcat(Projections4, [sum(abs.(eigenvecs[291:298, i]).^2) for i in 1:k])
        Projections5 = hcat(Projections5, [sum(abs.(eigenvecs[23:32, i]).^2) for i in 1:k])
        Projections6 = hcat(Projections6, [sum(abs.(eigenvecs[299:306, i]).^2) for i in 1:k])
        Projections7 = hcat(Projections7, [sum(abs.(eigenvecs[33:42, i]).^2) for i in 1:k])
        Projections8 = hcat(Projections8, [sum(abs.(eigenvecs[307:314, i]).^2) for i in 1:k])
    end

    h5open(output_filename, "w") do file
        file["eigenvalues"] = Energies

        file["proj_H_up"] = proj0_H_up
		file["proj_H_dn"] = proj0_H_dn
        file["proj1_dxy_up"] = proj1_dxy_up
        file["proj1_dxy_dn"] = proj1_dxy_dn
        file["proj1_dyz_up"] = proj1_dyz_up
        file["proj1_dyz_dn"] = proj1_dyz_dn
        file["proj1_dz2_up"] = proj1_dz2_up
        file["proj1_dz2_dn"] = proj1_dz2_dn
        file["proj1_dxz_up"] = proj1_dxz_up
        file["proj1_dxz_dn"] = proj1_dxz_dn
        file["proj1_x2y2_up"] = proj1_x2y2_up
        file["proj1_x2y2_dn"] = proj1_x2y2_dn
        file["proj2_sp1_up"] = proj2_sp1_up
        file["proj2_sp1_dn"] = proj2_sp1_dn
        file["proj2_sp2_up"] = proj2_sp2_up
        file["proj2_sp2_dn"] = proj2_sp2_dn
        file["proj2_sp3_up"] = proj2_sp3_up
        file["proj2_sp3_dn"] = proj2_sp3_dn
        file["proj2_sp4_up"] = proj2_sp4_up
        file["proj2_sp4_dn"] = proj2_sp4_dn

        file["proj1"] = Projections1
        file["proj2"] = Projections2
        file["proj3"] = Projections3
        file["proj4"] = Projections4
        file["proj5"] = Projections5
        file["proj6"] = Projections6
        file["proj7"] = Projections7
        file["proj8"] = Projections8
    end

end

#############################################################################
########## Processing the TB data for Collinear Spin with Hydrogen ##########
######################## on the Supercell 2 x 1 x 7 #########################
#############################################################################
function bandstructure_H(output_filename::String = joinpath(FILE_DATA_SOC_H_S,"BANDS_H.h5") )
    kx = vcat(zeros(100), range(0, 0.5, length = 50), zeros(100) .+ 0.5, range(0.5, 0, length = 111))  # Range of kx values
    ky = vcat(range(0.5, 0, length = 100), zeros(50),range(0, 0.5, length = 100),range(0.5, 0, length = 111)) # Range of ky values

    k = 1010
    N_tot = length(kx)

    # Compute the Hamiltonian for the given kx, ky pair

    eigvals, eigenvecs = eigen(bloch_hamiltonian(hr, SVector(kx[1], ky[1], 0.0)));

    # eigvals, sorted_indices = closest_elements(real(eigvals),  e_fermi, k)
    # eigenvecs = eigenvecs[:, sorted_indices]

    # Store the eigenvalues    
    Energies = sort(real.(eigvals))

    proj0_H_up = [abs.(eigenvecs[1,i]).^2 for i in 1:k]
	proj0_H_dn = [abs.(eigenvecs[2,i]).^2 for i in 1:k]

    Projections1 = [sum(abs.(eigenvecs[3:22, i]).^2) for i in 1:k]
    Projections2 = [sum(abs.(eigenvecs[563:578, i]).^2) for i in 1:k]
    Projections3 = [sum(abs.(eigenvecs[23:42, i]).^2) for i in 1:k]
    Projections4 = [sum(abs.(eigenvecs[579:594, i]).^2) for i in 1:k]
    Projections5 = [sum(abs.(eigenvecs[43:62, i]).^2) for i in 1:k]
    Projections6 = [sum(abs.(eigenvecs[595:610, i]).^2) for i in 1:k]
    Projections7 = [sum(abs.(eigenvecs[63:82, i]).^2) for i in 1:k]
    Projections8 = [sum(abs.(eigenvecs[611:626, i]).^2) for i in 1:k]
    for i in 2:length(kx)
        # Perform eigenvalue decomposition
        eigvals, eigenvecs = eigen(bloch_hamiltonian(hr, SVector(kx[i], ky[i], 0.0)));
        # eigvals, sorted_indices = closest_elements(real(eigvals),  e_fermi, k)
        # eigenvecs = eigenvecs[:, sorted_indices]

        # Store the eigenvalues
        Energies = hcat(Energies, sort(real.(eigvals)))

        proj0_H_up = hcat(proj0_H_up, [abs.(eigenvecs[1,i]).^2 for i in 1:k])
        proj0_H_dn = hcat(proj0_H_dn, [abs.(eigenvecs[2,i]).^2 for i in 1:k])
        Projections1 = hcat(Projections1, [sum(abs.(eigenvecs[3:22, i]).^2) for i in 1:k])
        Projections2 = hcat(Projections2, [sum(abs.(eigenvecs[563:578, i]).^2) for i in 1:k])
        Projections3 = hcat(Projections3, [sum(abs.(eigenvecs[23:42, i]).^2) for i in 1:k])
        Projections4 = hcat(Projections4, [sum(abs.(eigenvecs[579:594, i]).^2) for i in 1:k])
        Projections5 = hcat(Projections5, [sum(abs.(eigenvecs[43:62, i]).^2) for i in 1:k])
        Projections6 = hcat(Projections6, [sum(abs.(eigenvecs[595:610, i]).^2) for i in 1:k])
        Projections7 = hcat(Projections7, [sum(abs.(eigenvecs[63:82, i]).^2) for i in 1:k])
        Projections8 = hcat(Projections8, [sum(abs.(eigenvecs[611:626, i]).^2) for i in 1:k])
    end

    h5open(output_filename, "w") do file
        file["eigenvalues"] = Energies

        file["proj_H_up"] = proj0_H_up
		file["proj_H_dn"] = proj0_H_dn

        file["proj1"] = Projections1
        file["proj2"] = Projections2
        file["proj3"] = Projections3
        file["proj4"] = Projections4
        file["proj5"] = Projections5
        file["proj6"] = Projections6
        file["proj7"] = Projections7
        file["proj8"] = Projections8
    end

end

#############################################################################
######################### Bandstructure Plots ################################


#############################################################################
###################### Plots section: TB vs VASP ############################
#############################################################################
function TBvsVASP_nSOC_H()
    fermi_energy = 5.53063566
    
    fid = h5open(joinpath(FILE_DATA_nSOC_H,"BANDS_H.h5"), "r")

    Energies_up = read(fid, "eigenvalues_up")
    Energies_dn = read(fid, "eigenvalues_dn")
    Energies_up .= Energies_up .- fermi_energy
    Energies_dn .= Energies_dn .- fermi_energy
    # Initialize arrays to store reconstructed data
    proj1_dz2_up  = read(fid, "proj1_dxy_up")
    proj1_dz2_dn  = read(fid, "proj1_dxy_dn")
    proj1_dxz_up  = read(fid, "proj1_dyz_up")
    proj1_dxz_dn  = read(fid, "proj1_dyz_dn")
    proj1_dyz_up  = read(fid, "proj1_dz2_up")
    proj1_dyz_dn  = read(fid, "proj1_dz2_dn")
    proj1_x2y2_up = read(fid, "proj1_dxz_up")
    proj1_x2y2_dn = read(fid, "proj1_dxz_dn")
    proj1_dxy_up  = read(fid, "proj1_x2y2_up")
    proj1_dxy_dn  = read(fid, "proj1_x2y2_dn")
    
    proj1_up   = read(fid, "proj1_up")
    proj1_dn   = read(fid, "proj1_dn")
    proj2_up   = read(fid, "proj2_up")
    proj2_dn   = read(fid, "proj2_dn")
    proj3_up   = read(fid, "proj3_up")
    proj3_dn   = read(fid, "proj3_dn")
    proj4_up   = read(fid, "proj4_up")
    proj4_dn   = read(fid, "proj4_dn")



    close(fid)

    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    x_axis = range(0,1,size(Energies_up)[2])
    # Energies_1 = sort(enn, dims=1)
    LINE_COLOR = [:green, :red, :blue, :black]


    data_x = []
    data_y = []
    data_z = []
    k_x_tmp = 0.0
    for i in joinpath.(FILE_DATA_nSOC_H,readdir(FILE_DATA_nSOC_H)[endswith.(readdir(FILE_DATA_nSOC_H), "_LAYER1_UP.dat")])
        data_tmp = readdlm(i)[3:end,1:3]
        
        tmp = setdiff(1:size(data_tmp)[1], findall(data_tmp[:,1] .== "#" ) )
        
        data_x = vcat(data_x, data_tmp[:,1][tmp,:] .+ k_x_tmp)
        data_y = vcat(data_y, data_tmp[:,2][tmp,:])
        data_z = vcat(data_z, data_tmp[:,3][tmp,:])
        k_x_tmp += maximum( data_tmp[:,1][tmp,:])
    end
    data_x = collect((x for x in data_x))
    data_y = collect((x for x in data_y))
    data_z = collect((x for x in data_z))


    scatter!(p1, data_x./maximum(data_x), data_y, markersize=5.0*data_z, label=false, markercolor = :green, seriesalpha = 0.5)

    data_x = []
    data_y = []
    data_z = []
    k_x_tmp = 0.0
    for i in joinpath.(FILE_DATA_nSOC_H,readdir(FILE_DATA_nSOC_H)[endswith.(readdir(FILE_DATA_nSOC_H), "_LAYER1_DW.dat")])
        data_tmp = readdlm(i)[3:end,1:3]
        
        tmp = setdiff(1:size(data_tmp)[1], findall(data_tmp[:,1] .== "#" ) )
        
        data_x = vcat(data_x, data_tmp[:,1][tmp,:] .+ k_x_tmp)
        data_y = vcat(data_y, data_tmp[:,2][tmp,:])
        data_z = vcat(data_z, data_tmp[:,3][tmp,:])
        k_x_tmp += maximum( data_tmp[:,1][tmp,:])
    end
    data_x = collect((x for x in data_x))
    data_y = collect((x for x in data_y))
    data_z = collect((x for x in data_z))


    scatter!(p1, data_x./maximum(data_x), data_y, markersize=5.0*data_z, label=false, markercolor = :green, seriesalpha = 0.5)

    hline!(p1,[0],
        lw=LINE_WIDTH[3],
        lc=LINE_COLOR[3],
        ls=:dashdot,
        ylims=(-1.2,0.2),
        label=false
    )


    en_cut = -0.00
    hline!(p1,[en_cut],
        lw=LINE_WIDTH[3]+1,
        lc=:black,
        ls=:dashdot,
        ylims=(-1.2,0.2),
        label=false
    )


    vline!(p1,[0, 50/220, 100/220, 150/220, 220/220], lw=LINE_WIDTH[4],
    lc=:black,
    ls=:solid,label=false
    )
    plot!(p1, xticks = ([0, 100/440, 200/440, 300/440, 440/440],["Y", "\$ \\Gamma \$", "X", "M", "\$ \\Gamma \$"]),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12)

    plot!(p1,xlims=(0,1))
    
    for j in 1:size(Energies_dn)[1]
        plot!(p1, x_axis,Energies_up[j,:], label=false, linealpha = 0.3, linecolor = 1)
    end
    for j in 1:size(Energies_dn)[1]
        plot!(p1, x_axis,Energies_dn[j,:], label=false, linealpha = 0.3, linecolor = 1)
    end

    return p1

end

function TBvsVASP_SOC_H()
    fermi_energy = 5.5337
    k = 1:506
    p1 = plot()

    fid = h5open(joinpath(PATH_FILES_H,"BANDS_H.h5"), "r")

    Energies = read(fid, "eigenvalues")
    Energies .= Energies .- fermi_energy

    # Initialize arrays to store reconstructed data
    proj0_H_up = read(fid, "proj_H_up")
    proj0_H_dn = read(fid, "proj_H_dn")
    
    proj1   = read(fid, "proj1")
    proj2   = read(fid, "proj2")
    proj3   = read(fid, "proj3")
    proj4   = read(fid, "proj4")
    proj5   = read(fid, "proj5")
    proj6   = read(fid, "proj6")
    proj7   = read(fid, "proj7")
    proj8   = read(fid, "proj8")

    close(fid)

    length_k = size(Energies)[2]
    x_axis = range(0,1,length_k)

    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:green, :red, :blue, :black]
    for j in 1:length_k
        plot!(p1, x_axis, Energies[j,:], label=false, linealpha = 0.3, linecolor = 1)
    end

    # FILE_BANDS= raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/nCollinear_H_bridge/04_Bandstructure"

    data_x = []
    data_y = []
    data_z = []
    k_x_tmp = 0.0
    for i in joinpath.(PATH_FILES_H,readdir(PATH_FILES_H)[endswith.(readdir(PATH_FILES_H), "LAYER1.dat")])
        data_tmp = readdlm(i)[3:end,1:3]
        
        tmp = setdiff(1:size(data_tmp)[1], findall(data_tmp[:,1] .== "#" ) )
        
        data_x = vcat(data_x, data_tmp[:,1][tmp,:] .+ k_x_tmp)
        data_y = vcat(data_y, data_tmp[:,2][tmp,:])
        data_z = vcat(data_z, data_tmp[:,3][tmp,:])
        k_x_tmp += maximum( data_tmp[:,1][tmp,:])
    end

    data_x = collect((x for x in data_x))
    data_y = collect((x for x in data_y))
    data_z = collect((x for x in data_z))


    scatter!(p1, data_x./maximum(data_x), data_y, markersize=5.0 .* data_z, label=false, markercolor = :green, makeralpha = 0.5 .* data_z,
    markerstrokecolor = :black, markerstrokealpha = 0.0, markerstrokewidth = 0.0
    )


    hline!(p1,[0],
        lw=LINE_WIDTH[3],
        lc=LINE_COLOR[3],
        ls=:dashdot,
        ylims=(-1.2,0.2),
        label=false
    )


    en_cut = -0.00
    hline!(p1,[en_cut],
        lw=LINE_WIDTH[3]+1,
        lc=:black,
        ls=:dashdot,
        ylims=(-1.2,0.2),
        label=false
    )


    vline!(p1,[0, 50/220, 100/220, 150/220, 220/220], lw=LINE_WIDTH[4],
    lc=:black,
    ls=:solid,label=false
    )
    plot!(p1, xticks = ([0, 100/440, 200/440, 300/440, 440/440],["Y", "\$ \\Gamma \$", "X", "M", "\$ \\Gamma \$"]),
        ytickfontsize  = 10,
        yticks = -1.2:0.2:0.2,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12)

    plot!(p1,xlims=(0,1))
    
    # savefig("/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell/Figures/Bandstructure/New/TB_vs_VASP_surf_H.png")


    return p1

end

function TBvsVASP_S217_H()
    fermi_energy = 5.7468
    
    fid = h5open(joinpath(PATH_FILES_H_S,"BANDS_H.h5"), "r")

    Energies = read(fid, "eigenvalues")
    Energies .= Energies .- fermi_energy

    # Initialize arrays to store reconstructed data
    proj0_H_up = read(fid, "proj_H_up")
    proj0_H_dn = read(fid, "proj_H_dn")
    
    proj1   = read(fid, "proj1")
    proj2   = read(fid, "proj2")
    proj3   = read(fid, "proj3")
    proj4   = read(fid, "proj4")
    proj5   = read(fid, "proj5")
    proj6   = read(fid, "proj6")
    proj7   = read(fid, "proj7")
    proj8   = read(fid, "proj8")

    close(fid)

    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    x_axis = range(0,1,size(Energies)[2])
    # Energies_1 = sort(enn, dims=1)
    LINE_COLOR = [:green, :red, :blue, :black]


    data_x = []
    data_y = []
    data_z = []
    k_x_tmp = 0.0
    for i in joinpath.(PATH_FILES_H_S,readdir(PATH_FILES_H_S)[endswith.(readdir(PATH_FILES_H_S), "_LAYER1.dat")])
        data_tmp = readdlm(i)[3:end,1:3]
        
        tmp = setdiff(1:size(data_tmp)[1], findall(data_tmp[:,1] .== "#" ) )
        
        data_x = vcat(data_x, data_tmp[:,1][tmp,:] .+ k_x_tmp)
        data_y = vcat(data_y, data_tmp[:,2][tmp,:])
        data_z = vcat(data_z, data_tmp[:,3][tmp,:])
        k_x_tmp += maximum( data_tmp[:,1][tmp,:])
    end
    data_x = collect((x for x in data_x))
    data_y = collect((x for x in data_y))
    data_z = collect((x for x in data_z))


    scatter!(p1, data_x./maximum(data_x), data_y, markersize=5.0*data_z, label=false, markercolor = :green, seriesalpha = 0.5)


    hline!(p1,[0],
        lw=LINE_WIDTH[3],
        lc=LINE_COLOR[3],
        ls=:dashdot,
        ylims=(-1.2,0.2),
        label=false
    )


    en_cut = -0.00
    hline!(p1,[en_cut],
        lw=LINE_WIDTH[3]+1,
        lc=:black,
        ls=:dashdot,
        ylims=(-1.2,0.2),
        label=false
    )


    vline!(p1,[0, 100/361, 150/361, 250/361, 361/361], lw=LINE_WIDTH[4],
    lc=:black,
    ls=:solid,label=false
    )
    plot!(p1, xticks = ([0, 100/361, 150/361, 250/361, 361/361],["Y", "\$ \\Gamma \$", "X", "M", "\$ \\Gamma \$"]),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12)

    plot!(p1,xlims=(0,1))
    
    for j in 1:size(Energies)[1]
        plot!(p1, x_axis,Energies[j,:], label=false, linealpha = 0.3, linecolor = 1)
    end

    return p1
end

#############################################################################
####################### Plots section: TB PBANDS ############################
#############################################################################

function bandstructure_nSOC()
    fermi_energy = 5.58199864
    k = 1:506
    p1 = plot()

    fid = h5open(joinpath(FILE_DATA_nSOC,"BANDS.h5"), "r")

    Energies_up = read(fid, "eigenvalues_up")
    Energies_dn = read(fid, "eigenvalues_dn")
    Energies_up .= Energies_up .- fermi_energy
    Energies_dn .= Energies_dn .- fermi_energy

    # Initialize arrays to store reconstructed data
    proj_up   = read(fid, "proj1_up") .+ read(fid, "proj2_up") .+ read(fid, "proj3_up") .+ read(fid, "proj4_up") 
    proj_dn   = read(fid, "proj1_dn") .+ read(fid, "proj2_dn") .+ read(fid, "proj3_dn") .+ read(fid, "proj4_dn") 


    close(fid)

    length_k = size(Energies_up)[2]
    x_axis = range(0,1,length_k)
    
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:black, :blue, :red, :green]



    plot!( p1,
            ylabel       = "Energy " * L"- E_F" * " (eV) ",
            title        = "Electronic Band Structure",
            legend       = false,
            framestyle   = :box,
            ygrid        = false,
            ylim         = (-1.2, 0.2),
            xlim        = (0.0, 1.0)
        )
    plot!(p1,
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid          = false
    )
    vline!(p1,
        [0, 50/220, 100/220, 150/220, 220/220],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )
    hline!(p1,
        [0.0],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )
    plot!(p1, 
        xticks = ([0, 100/440, 200/440, 300/440, 440/440],["Y", "\$ \\Gamma \$", "X", "M", "\$ \\Gamma \$"]),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12
    )

    for j in 1:size(Energies_up)[1]
        plot!(
            p1,
            x_axis,
            Energies_up[j,:],
            label     = false,
            linealpha = 0.3,
            ls        = :solid,
            linecolor = LINE_COLOR[1]
        )
        plot!(
            p1,
            x_axis,
            Energies_dn[j,:],
            label     = false,
            linealpha = 0.3,
            ls        = :solid,
            linecolor = LINE_COLOR[1]
        )
    end
    Markersize = [4.0 5.0]
    for j in 1:size(Energies_up)[1]
        scatter!(
            p1,
            x_axis,
            Energies_up[j,:],
            label             = false,
            markersize        = Markersize[1] .* proj_up[j,:],
            markerstrokewidth = 0.2,
            markercolor       = 1,
            markeralpha       = Markersize[1] .* proj_up[j,:]
        )
        scatter!(
            p1,
            x_axis,
            Energies_dn[j,:],
            label             = false,
            markersize        = Markersize[1] .* proj_dn[j,:],
            markerstrokewidth = 0.2,
            markercolor       = 1,
            markeralpha       = Markersize[1] .* proj_dn[j,:]
        )
    end
    # savefig("/media/r_floren/FernandoM/VASP/NbP/Supercell/Figures/Bandstructure/New/Slab_H/TB/Bands_surf.png")
    return p1

end


function bandstructure_nSOC_H()
    fermi_energy = 5.53063566
    k = 1:506
    p1 = plot()

    fid = h5open(joinpath(FILE_DATA_nSOC_H,"BANDS_H.h5"), "r")

    Energies_up = read(fid, "eigenvalues_up")
    Energies_dn = read(fid, "eigenvalues_dn")
    Energies_up .= Energies_up .- fermi_energy
    Energies_dn .= Energies_dn .- fermi_energy

    # Initialize arrays to store reconstructed data
    proj0_H_up = read(fid, "proj_H_up")
    proj0_H_dn = read(fid, "proj_H_dn")
    
    proj_up   = read(fid, "proj1_up") .+ read(fid, "proj2_up") .+ read(fid, "proj3_up") .+ read(fid, "proj4_up") 
    proj_dn   = read(fid, "proj1_dn") .+ read(fid, "proj2_dn") .+ read(fid, "proj3_dn") .+ read(fid, "proj4_dn") 


    close(fid)

    length_k = size(Energies_up)[2]
    x_axis = range(0,1,length_k)
    
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:black, :blue, :red, :green]



    plot!( p1,
            ylabel       = "Energy " * L"- E_F" * " (eV) ",
            title        = "Electronic Band Structure",
            legend       = false,
            framestyle   = :box,
            ygrid        = false,
            ylim         = (-1.2, 0.2),
            xlim        = (0.0, 1.0)
        )
    plot!(p1,
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid          = false
    )
    vline!(p1,
        [0, 50/220, 100/220, 150/220, 220/220],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )
    hline!(p1,
        [0.0],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )
    plot!(p1, 
        xticks = ([0, 100/440, 200/440, 300/440, 440/440],["Y", "\$ \\Gamma \$", "X", "M", "\$ \\Gamma \$"]),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12
    )

    for j in 1:size(Energies_up)[1]
        plot!(
            p1,
            x_axis,
            Energies_up[j,:],
            label     = false,
            linealpha = 0.3,
            ls        = :solid,
            linecolor = LINE_COLOR[1]
        )
        plot!(
            p1,
            x_axis,
            Energies_dn[j,:],
            label     = false,
            linealpha = 0.3,
            ls        = :solid,
            linecolor = LINE_COLOR[1]
        )
    end
    Markersize = [4.0 5.0]
    for j in 1:size(Energies_up)[1]
        scatter!(
            p1,
            x_axis,
            Energies_up[j,:],
            label             = false,
            markersize        = Markersize[1] .* proj_up[j,:],
            markerstrokewidth = 0.2,
            markercolor       = 1,
            markeralpha       = Markersize[1] .* proj_up[j,:]
        )
        scatter!(
            p1,
            x_axis,
            Energies_up[j,:],
            label             = false,
            markersize        = 2.0, #Markersize[2] .* (proj0_H[j,:]),
            markerstrokewidth = 0.2,
            markercolor       = LINE_COLOR[3],
            markeralpha       = Markersize[2] .* (proj0_H_up[j,:])
        )
        scatter!(
            p1,
            x_axis,
            Energies_dn[j,:],
            label             = false,
            markersize        = Markersize[1] .* proj_dn[j,:],
            markerstrokewidth = 0.2,
            markercolor       = 1,
            markeralpha       = Markersize[1] .* proj_dn[j,:]
        )
        scatter!(
            p1,
            x_axis,
            Energies_dn[j,:],
            label             = false,
            markersize        = 2.0, #Markersize[2] .* (proj0_H[j,:]),
            markerstrokewidth = 0.2,
            markercolor       = LINE_COLOR[3],
            markeralpha       = Markersize[2] .* (proj0_H_dn[j,:])
        )
    end
    # savefig("/media/r_floren/FernandoM/VASP/NbP/Supercell/Figures/Bandstructure/New/Slab_H/TB/Bands_surf.png")
    return p1

end

function bandstructure_SOC()
    fermi_energy = 5.62997313
    k = 1:506
    p1 = plot()

    fid = h5open(joinpath(FILE_DATA_SOC,"BANDS.h5"), "r")

    Energies  = read(fid, "eigenvalues")
    Energies .= Energies .- fermi_energy

    # Initialize arrays to store reconstructed data
    proj   = read(fid, "proj1") .+ read(fid, "proj2") .+ read(fid, "proj3") .+ read(fid, "proj4") 

    close(fid)

    length_k = size(Energies)[2]
    x_axis = range(0,1,length_k)
    
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:black, :blue, :red, :green]



    plot!( p1,
            ylabel       = "Energy " * L"- E_F" * " (eV) ",
            title        = "Electronic Band Structure",
            legend       = false,
            framestyle   = :box,
            ygrid        = false,
            ylim         = (-2.2, 1.2),
            xlim        = (0.0, 1.0)
        )
    plot!(p1,
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid          = false
    )
    hline!(p1,
        [0.0],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )


    for j in 1:size(Energies)[1]
        plot!(
            p1,
            x_axis,
            Energies[j,:],
            label     = false,
            linealpha = 0.3,
            ls        = :solid,
            linecolor = LINE_COLOR[1]
        )
    end
    Markersize = [4.0 7.0]
    for j in 1:size(Energies)[1]
        scatter!(
            p1,
            x_axis,
            Energies[j,:],
            label             = false,
            markersize        = Markersize[1] .* proj[j,:],
            markerstrokewidth = 0.2,
            markercolor       = 1,
            markeralpha       = Markersize[1] .* proj[j,:]
        )
    end
    # savefig("/media/r_floren/FernandoM/VASP/NbP/Supercell/Figures/Bandstructure/New/Slab_H/TB/Bands_surf.png")
    return p1

end

function bandstructure_SOC_tst()
    fermi_energy = 5.58485570
    k = 1:504
    p1 = plot()

    fid = h5open(joinpath(@__DIR__,"BANDS_gpu.h5"), "r")

    Energies  = read(fid, "eigenvalues")
    Energies .= Energies .- fermi_energy

    # Initialize arrays to store reconstructed data
    proj   = read(fid, "proj1") .+ read(fid, "proj2") .+ read(fid, "proj3") .+ read(fid, "proj4") 
    # proj   = read(fid, "proj3")

    close(fid)

    length_k = size(Energies)[2]
    x_axis = range(0,1,length_k)
    
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:black, :blue, :red, :green]



    plot!( p1,
            ylabel       = "Energy " * L"- E_F" * " (eV) ",
            title        = "Electronic Band Structure",
            legend       = false,
            framestyle   = :box,
            ygrid        = false,
            ylim         = (-1.2, 0.2),
            xlim        = (0.0, 1.0)
        )
    plot!(p1,
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid          = false
    )
    
    hline!(p1,
        [0.0],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )

    for j in 1:size(Energies)[1]
        plot!(
            p1,
            x_axis,
            Energies[j,:],
            label     = false,
            linealpha = 0.3,
            ls        = :solid,
            linecolor = LINE_COLOR[1]
        )
    end
    Markersize = [7.0 10.0]
    for j in 1:size(Energies)[1]
        scatter!(
            p1,
            x_axis,
            Energies[j,:],
            label             = false,
            markersize        = Markersize[2] .* proj[j,:],
            markerstrokewidth = 0.2,
            markercolor       = 1,
            markeralpha       = Markersize[1] .* proj[j,:]
        )
    end
    # savefig("/media/r_floren/FernandoM/VASP/NbP/Supercell/Figures/Bandstructure/New/Slab_H/TB/Bands_surf.png")
    return p1

end

function bandstructure_SOC_H()
    fermi_energy = 5.53348380
    k = 1:506
    p1 = plot()

    fid = h5open(joinpath(FILE_DATA_SOC_H2,"BANDS_H2.h5"), "r")

    Energies  = read(fid, "eigenvalues")
    Energies .= Energies .- fermi_energy

    # Initialize arrays to store reconstructed data
    proj0_H = read(fid, "proj_H_up") .+ read(fid, "proj_H_dn")
    
    proj   = read(fid, "proj1") .+ read(fid, "proj2") .+ read(fid, "proj3") .+ read(fid, "proj4") 

    close(fid)

    length_k = size(Energies)[2]
    x_axis = range(0,1,length_k)
    
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:black, :blue, :red, :green]



    plot!( p1,
            ylabel       = "Energy " * L"- E_F" * " (eV) ",
            title        = "Electronic Band Structure",
            legend       = false,
            framestyle   = :box,
            ygrid        = false,
            ylim         = (-1.2, 0.2),
            xlim        = (0.0, 1.0)
        )
    plot!(p1,
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid          = false
    )
    vline!(p1,
        [0, 50/220, 100/220, 150/220, 220/220],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )
    hline!(p1,
        [0.0],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )
    plot!(p1, 
        xticks = ([0, 100/440, 200/440, 300/440, 440/440],["Y", "\$ \\Gamma \$", "X", "M", "\$ \\Gamma \$"]),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12
    )

    for j in 1:size(Energies)[1]
        plot!(
            p1,
            x_axis,
            Energies[j,:],
            label     = false,
            linealpha = 0.3,
            ls        = :solid,
            linecolor = LINE_COLOR[1]
        )
    end
    Markersize = [4.0 5.0]
    for j in 1:size(Energies)[1]
        scatter!(
            p1,
            x_axis,
            Energies[j,:],
            label             = false,
            markersize        = Markersize[1] .* proj[j,:],
            markerstrokewidth = 0.2,
            markercolor       = 1,
            markeralpha       = Markersize[1] .* proj[j,:]
        )
        scatter!(
            p1,
            x_axis,
            Energies[j,:],
            label             = false,
            markersize        = 2.0, #Markersize[2] .* (proj0_H[j,:]),
            markerstrokewidth = 0.2,
            markercolor       = LINE_COLOR[3],
            markeralpha       = Markersize[1] .* (proj0_H[j,:])
        )
    end
    # savefig("/media/r_floren/FernandoM/VASP/NbP/Supercell/Figures/Bandstructure/New/Slab_H/TB/Bands_surf.png")
    return p1

end

function bandstructure_SOC_H_tst()
    fermi_energy = 5.53348380
    k = 1:506
    p1 = plot()

    fid = h5open(joinpath(@__DIR__,"BANDS_H2.h5"), "r")

    Energies  = read(fid, "eigenvalues")
    Energies .= Energies .- fermi_energy

    # Initialize arrays to store reconstructed data
    proj0_H = read(fid, "proj_H_up") .+ read(fid, "proj_H_dn")
    
    proj   = read(fid, "proj1") .+ read(fid, "proj2") .+ read(fid, "proj3") .+ read(fid, "proj4") 

    close(fid)

    length_k = size(Energies)[2]
    x_axis = range(0,1,length_k)
    
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:black, :blue, :red, :green]



    plot!( p1,
            ylabel       = "Energy " * L"- E_F" * " (eV) ",
            title        = "Electronic Band Structure",
            legend       = false,
            framestyle   = :box,
            ygrid        = false,
            ylim         = (-1.2, 0.2),
            xlim        = (0.0, 1.0)
        )
    plot!(p1,
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid          = false
    )
    hline!(p1,
        [0.0],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )


    for j in 1:size(Energies)[1]
        plot!(
            p1,
            x_axis,
            Energies[j,:],
            label     = false,
            linealpha = 0.3,
            ls        = :solid,
            linecolor = LINE_COLOR[1]
        )
    end
    Markersize = [7.0 10.0]
    for j in 1:size(Energies)[1]
        scatter!(
            p1,
            x_axis,
            Energies[j,:],
            label             = false,
            markersize        = Markersize[2] .* proj[j,:],
            markerstrokewidth = 0.2,
            markercolor       = 1,
            markeralpha       = Markersize[1] .* proj[j,:]
        )
        scatter!(
            p1,
            x_axis,
            Energies[j,:],
            label             = false,
            markersize        = 2.0, #Markersize[2] .* (proj0_H[j,:]),
            markerstrokewidth = 0.2,
            markercolor       = LINE_COLOR[3],
            markeralpha       = Markersize[1] .* (proj0_H[j,:])
        )
    end
    # savefig("/media/r_floren/FernandoM/VASP/NbP/Supercell/Figures/Bandstructure/New/Slab_H/TB/Bands_surf.png")
    return p1

end

function bandstructure_SOC_H_S()
    fermi_energy = 5.74715706
    p1 = plot()

    fid = h5open(joinpath(FILE_DATA_SOC_H_S,"BANDS_H.h5"), "r")

    Energies  = read(fid, "eigenvalues")
    Energies .= Energies .- fermi_energy

    # Initialize arrays to store reconstructed data
    proj0_H = read(fid, "proj_H_up") .+ read(fid, "proj_H_dn")
    
    proj   = read(fid, "proj1") .+ read(fid, "proj2") .+ read(fid, "proj3") .+ read(fid, "proj4") 

    close(fid)

    length_k = size(Energies)[2]
    x_axis = range(0,1,length_k)
    
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:black, :blue, :red, :green]

    plot!( p1,
            ylabel       = "Energy " * L"- E_F" * " (eV) ",
            title        = "Electronic Band Structure",
            legend       = false,
            framestyle   = :box,
            ygrid        = false,
            ylim         = (-1.2, 0.2),
            xlim        = (0.0, 1.0)
        )
    plot!(p1,
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid          = false
    )
    vline!(p1,
        [0, 50/181, 75/181, 125/181, 181/181],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )
    hline!(p1,
        [0.0],
        lw    = LINE_WIDTH[4],
        lc    = :black,
        ls    = :solid,
        label = false
    )
    plot!(p1, 
        xticks = ([0, 50/181, 75/181, 125/181, 181/181],["Y", "\$ \\Gamma \$", "X", "M", "\$ \\Gamma \$"]),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12
    )

    for j in 1:size(Energies)[1]
        plot!(
            p1,
            x_axis,
            Energies[j,:],
            label     = false,
            linealpha = 0.3,
            ls        = :solid,
            linecolor = LINE_COLOR[1]
        )
    end
    Markersize = [4.0 5.0]
    for j in 1:size(Energies)[1]
        scatter!(
            p1,
            x_axis,
            Energies[j,:],
            label             = false,
            markersize        = Markersize[1] .* proj[j,:],
            markerstrokewidth = 0.2,
            markercolor       = 1,
            markeralpha       = Markersize[1] .* proj[j,:]
        )
        scatter!(
            p1,
            x_axis,
            Energies[j,:],
            label             = false,
            markersize        = 2.0, #Markersize[2] .* (proj0_H[j,:]),
            markerstrokewidth = 0.2,
            markercolor       = LINE_COLOR[3],
            markeralpha       = Markersize[2] .* (proj0_H[j,:])
        )
    end
    # savefig("/media/r_floren/FernandoM/VASP/NbP/Supercell/Figures/Bandstructure/New/Slab_H/TB/Bands_surf.png")
    return p1

end


#############################################################################
####################### Hydrogen Bands section ############################
#############################################################################


"Read a seedname_hr.dat file and return an `HRData` object."
function load_hr_H(path::String)
    open(path, "r") do io
        nw_H = 2
        readline(io)                            # comment / timestamp line
        nw      = parse(Int, readline(io))      # N_R
        nr      = parse(Int, readline(io))      # N_R
        # read nr degeneracy integers — the file prints 15 per line
        w = Int[];  while length(w) < nr
            append!(w, parse.(Int, split(strip(readline(io)))))
        end
        # allocate container for R‑vectors and matrix elements
        Rvec = SVector(0,0,1)
        Rvec = vcat(Rvec',Rvec')
        for i in 1:nr-2
            Rvec = vcat(Rvec, SVector(0,0,1)')
        end
        Rvec = [Rvec[i,:] for i in 1:size(Rvec,1)]
        H_R   = Array{ComplexF64}(undef, nw_H, nw_H, nr)

        # every line:  Rx Ry Rz   i   j   Re   Im
        # indices i,j are 1‑based in the file → no shift needed
        tmp = 1
        allowed = Set([(1,2), (1,1), (2,1), (2,2)])
        for _ in 1:(nr * nw * nw)
            vals = split(readline(io))
            pair = (parse(Int, vals[4]), parse(Int, vals[5]))
            if pair in allowed
                Rx,Ry,Rz = parse.(Int, vals[1:3])
                ii, jj   = parse.(Int, vals[4:5])
                Re, Im   = parse.(Float64, vals[6]), parse(Float64, vals[7])
                rindex   = findfirst(r -> r == SVector(Rx,Ry,Rz), Rvec)
                if rindex === nothing                       # first time we meet this R
                    rindex = tmp                            # next empty slot
                    Rvec[rindex] = SVector(Rx,Ry,Rz)
                    tmp += 1
                end
                H_R[ii, jj, rindex] = complex(Re, Im)
            end
            
        end
        return HRData(nw_H, nr, Rvec, w, H_R)
    end
end

function SOC_H_bands()

    hr = load_hr_H(joinpath(FILE_DATA_SOC_H,"wannier90_hr.dat"))
    kx = range(-0.5, 0.5, length = 200) # Range of kx values
    ky = range(-0.5, 0.5, length = 200) # Range of kx values

    Energies = zeros(Float64, length(kx), length(ky), 2)



    for (idx,kxi) in enumerate(kx), (idy,kyi) in enumerate(ky)
        # Perform eigenvalue decomposition
        Energies[idx, idy, :] = eigen(bloch_hamiltonian(hr, SVector(kxi, kyi, 0.0))).values;

    end


    # plot(kx,ky,Energies[:,:,1],st=:surface,camera=(-50,30), alpha=0.5,
    #     xlabel = "kx", ylabel = "ky", zlabel = "Energy (eV)",
    #     title = "Hydrogen SOC Bands",
    #     legend = false, framestyle = :box
    # )
    # plot!(kx,ky,Energies[:,:,2],st=:surface,camera=(-30,30),alpha=0.5)

    return kx, ky, Energies

end