using LaTeXStrings
using Plots
using Distributions
using DelimitedFiles

const FILE_PROCAR_nSOC    = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/Collinear/05_V2W/out_v2w_dn_20250617_223128/PROCAR"
const FILE_PROCAR_nSOC_H  = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/Collinear_H/05_V2W/out_v2w_up_20250611_000304/PROCAR"
const FILE_PROCAR_SOC     = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/nCollinear/05_V2W/out_v2w_gpu_20250609_035612/PROCAR"
const FILE_PROCAR_SOC_H   = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/nCollinear_H/08_VASP_PP/DOS/out_nscf_20251117_114104/PROCAR"
const FILE_PROCAR_SOC_H_S = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell_217/nCollinear_H/05_V2W/out_v2w_20250607_000616/PROCAR"


# Fermi Dirac Distribution
function FD_dist(ϵ::Float64, T::Float64, E_fermi::Float64)
    return 1/(exp((ϵ - E_fermi)/T) + 1)
end

## Function to plot calculate the DOS and PDOS
"""
    PROCAREData
Container for the contents of a `seedname_hr.dat` file.
Fields
-------
  num_k      :: Int                      # number of k-points sampled in the BZ
  num_bnds   :: Int                      # number of bands
  num_ions   :: Int                      # number of atoms in the unit cell
  Energy     :: Matrix{Float64}          # Corresponding energi at each k-point with size (num_k,num_bnds)
  Projected  :: Matrix{Float64}          # size (num_k,num_bnds,num_ions)
  k_w        :: Vector{Float64}          # k-point weights, size (num_k,)
"""
struct PROCARDataDOS
    num_k       :: Int
    num_bnds    :: Int
    num_ions    :: Int
    Energy      :: Matrix{Float64}
    k_w         :: Vector{Float64}
end

struct PROCARDataPDOS
    num_k       :: Int
    num_bnds    :: Int
    num_ions    :: Int
    Energy      :: Matrix{Float64}
    Projected   :: Matrix{Float64}
    k_w         :: Vector{Float64}
end

struct PROCARDataPDOSH
    num_k       :: Int
    num_bnds    :: Int
    num_ions    :: Int
    Energy      :: Matrix{Float64}
    Projected   :: Matrix{Float64}
    ProjectedH  :: Matrix{Float64}
    k_w         :: Vector{Float64}
end

struct PROCARData_SOC_MAG
    num_k       :: Int
    num_bnds    :: Int
    num_ions    :: Int
    Energy      :: Matrix{Float64}
    Projected   :: Matrix{Float64}
    M_x         :: Matrix{Float64}  # Projected in the x direction
    M_y         :: Matrix{Float64}  # Projected in the y direction
    M_z         :: Matrix{Float64}  # Projected in the z direction
    k_w         :: Vector{Float64}
end
struct PROCARData_SOC_MAG_H
    num_k       :: Int
    num_bnds    :: Int
    num_ions    :: Int
    Energy      :: Matrix{Float64}
    Projected   :: Matrix{Float64}
    ProjectedH  :: Matrix{Float64}
    M_x         :: Matrix{Float64}  # Projected in the x direction
    M_y         :: Matrix{Float64}  # Projected in the y direction
    M_z         :: Matrix{Float64}  # Projected in the z direction
    k_w         :: Vector{Float64}
end

"Read a PROCAR file and return an `PROCARData` object."
function load_PROCAR_DOS_nSOC(path::String)
    num_kpoints = []; num_bnds = []; num_ions = []; Energies = []; k_weight = [];
    open(path, "r") do io
        readline(io)                            # comment / timestamp line
        tmp = split(readline(io))
        num_kpoints = parse(Int,tmp[4])
        num_bnds    = parse(Int,tmp[8])
        num_ions    = parse(Int,tmp[12])
        Energies    = zeros(2*num_kpoints, num_bnds)
        k_weight    = zeros(2*num_kpoints)
        readline(io)
        for i in 1:num_kpoints
            k_weight[i] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = split(readline(io))
                Energies[i, j] = parse(Float64, tmp[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                for i in 1:2+(1+num_ions)-1
                    readline(io)  # Skip the rest of the lines for this k-point
                end
            end
            readline(io)
        end
        readline(io)  # Read the last line (comment line)
        for i in 1:num_kpoints
            k_weight[i+num_kpoints] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = split(readline(io))
                Energies[i + num_kpoints, j] = parse(Float64, tmp[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                for i in 1:2+(1+num_ions)-1
                    readline(io)  # Skip the rest of the lines for this k-point
                end
            end
            readline(io)
        end  
    end

    return PROCARDataDOS(2*num_kpoints, num_bnds, num_ions, Energies, k_weight)
end

function load_PROCAR_PDOS_nSOC(path::String)
    num_kpoints = []; num_bnds = []; num_ions = []; Energies = []; PDOS = []; k_weight = [];
    open(path, "r") do io
        readline(io)                            # comment / timestamp line
        tmp = split(readline(io))
        num_kpoints = parse(Int,tmp[4])
        num_bnds    = parse(Int,tmp[8])
        num_ions    = parse(Int,tmp[12])
        Energies    = zeros(2*num_kpoints, num_bnds)
        PDOS        = zeros(2*num_kpoints, num_bnds)  # 3 for s, p, d orbitals and (1+num_ions)*4 for each ion's contributions
        k_weight    = zeros(2*num_kpoints)
        readline(io)
        for i in 1:num_kpoints
            k_weight[i] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = split(readline(io))
                Energies[i, j] = parse(Float64, tmp[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                for i in 1:num_ions
                    tmp = readline(io)
                    if parse(Int64, split(tmp)[1]) == 1 || parse(Int64, split(tmp)[1]) == 2 || parse(Int64, split(tmp)[1]) == 29 || parse(Int64, split(tmp)[1]) == 30 
                        PDOS[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                    end
                end
                readline(io)
                readline(io)
            end
            readline(io)
        end
        readline(io)  # Read the last line (comment line)

        for i in 1:num_kpoints
            k_weight[i+num_kpoints] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = split(readline(io))
                Energies[i + num_kpoints, j] = parse(Float64, tmp[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                for i in 1:num_ions
                    tmp = readline(io)
                    if parse(Int64, split(tmp)[1]) == 1 || parse(Int64, split(tmp)[1]) == 2 || parse(Int64, split(tmp)[1]) == 29 || parse(Int64, split(tmp)[1]) == 30 
                        PDOS[i + num_kpoints, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                    end
                end
                readline(io)  # Skip the rest of the lines for this k-point
                readline(io)
            end
            readline(io)
        end  
    end

    return PROCARDataPDOS(2*num_kpoints, num_bnds, num_ions, Energies, PDOS, k_weight)
end

function load_PROCAR_PDOS_nSOC_H(path::String)
    num_kpoints = []; num_bnds = []; num_ions = []; Energies = [];
    PDOS = []; PDOS_H = []; k_weight = [];
    open(path, "r") do io
        readline(io)                            # comment / timestamp line
        tmp = split(readline(io))
        num_kpoints = parse(Int,tmp[4])
        num_bnds    = parse(Int,tmp[8])
        num_ions    = parse(Int,tmp[12])
        Energies    = zeros(2*num_kpoints, num_bnds)
        PDOS        = zeros(2*num_kpoints, num_bnds)  # 3 for s, p, d orbitals and (1+num_ions)*4 for each ion's contributions
        PDOS_H      = zeros(2*num_kpoints, num_bnds)  # 3 for s, p, d orbitals and (1+num_ions)*4 for each ion's contributions
        k_weight    = zeros(2*num_kpoints)
        readline(io)
        for i in 1:num_kpoints
            k_weight[i] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = split(readline(io))
                Energies[i, j] = parse(Float64, tmp[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                for _ in 1:num_ions
                    tmp = readline(io)
                    # PDOS for the Hydrogen atom
                    if parse(Int64, split(tmp)[1]) == 1 
                        PDOS_H[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                    end
                    # PDOS for the surfacae atoms
                    if parse(Int64, split(tmp)[1]) == 2 || parse(Int64, split(tmp)[1]) == 3 || parse(Int64, split(tmp)[1]) == 4 || parse(Int64, split(tmp)[1]) == 5 || parse(Int64, split(tmp)[1]) == 30 || parse(Int64, split(tmp)[1]) == 31 || parse(Int64, split(tmp)[1]) == 32 || parse(Int64, split(tmp)[1]) == 33
                        PDOS[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                    end
                end
                readline(io)
                readline(io)
            end
            readline(io)
        end
        readline(io)  # Read the last line (comment line)

        for i in 1:num_kpoints
            k_weight[i+num_kpoints] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = split(readline(io))
                Energies[i + num_kpoints, j] = parse(Float64, tmp[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                for _ in 1:num_ions
                    tmp = readline(io)
                    # PDOS for the Hydrogen atom
                    if parse(Int64, split(tmp)[1]) == 1 
                        PDOS_H[i + num_kpoints, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                    end
                    # PDOS for the surfacae atoms
                    if parse(Int64, split(tmp)[1]) == 2 || parse(Int64, split(tmp)[1]) == 3 || parse(Int64, split(tmp)[1]) == 4 || parse(Int64, split(tmp)[1]) == 5 || parse(Int64, split(tmp)[1]) == 30 || parse(Int64, split(tmp)[1]) == 31 || parse(Int64, split(tmp)[1]) == 32 || parse(Int64, split(tmp)[1]) == 33
                        PDOS[i + num_kpoints, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                    end
                end
                readline(io)  # Skip the rest of the lines for this k-point
                readline(io)
            end
            readline(io)
        end  
    end

    return PROCARDataPDOSH(2*num_kpoints, num_bnds, num_ions, Energies, PDOS, PDOS_H, k_weight)
end


" This function computes the DOS using the provided data obtained from VASPKIT."

function DOS_VASPKIT()
    data_tmp = readdlm(joinpath(PATH_FILES_H, "PDOS_A1_SOC.dat"))
    e_fermi     = 5.533484
    data_energy = data_tmp[2:end,1]
    data_s      = data_tmp[2:end,2]
    data_py     = data_tmp[2:end,3]
    data_pz     = data_tmp[2:end,4]
    data_px     = data_tmp[2:end,5]
    data_tot    = data_tmp[2:end,11]

    p = findfirst(x -> x > e_fermi, data_energy)
    println(sum(data_s[1:p-1]) / sum(data_tot[1:p-1]))

    ##############################################################################
    # 2. DOS-style smearing
    ##############################################################################
    # --- user parameters ---
    σ        = 0.09  # smearing width
    Emin     = minimum(data_energy)
    Emax     = maximum(data_energy)
    nbins    = 2000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    pdos_s   = zeros(Float64, nbins)
    pdos_py  = zeros(Float64, nbins)
    pdos_pz  = zeros(Float64, nbins)
    pdos_px  = zeros(Float64, nbins)
    for idx in 1:length(data_energy)
        Ei   = data_energy[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            gval       = gnorm * pdf(gauss, Ei - ebins[j])
            pdos_s[j]  += gval * data_s[idx]
            pdos_py[j] += gval * data_py[idx]
            pdos_pz[j] += gval * data_pz[idx]
            pdos_px[j] += gval * data_px[idx]
        end
    end
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:green, :blue, :black]
    plot!(p1, ebins, pdos_s,  label=L"H_{s}", lc = 1, ls = :solid)
    plot!(p1, ebins, pdos_px, label=L"H_{px}", lc = 2, ls = :dash)
    plot!(p1, ebins, pdos_py, label=L"H_{py}", lc = 3, ls = :dash)
    plot!(p1, ebins, pdos_pz, label=L"H_{pz}", lc = 4, ls = :dash)
    plot!(p1,
        xlims=(-10,10),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12
    )
    return p1
end

####################################################################
######################### nSOC #####################################
####################################################################
" This function computes the DOS using the provided data obtained directly from the PROCAR file."
function DOS_PROCAR_nSOC()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.58199864
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_DOS_nSOC( FILE_PROCAR_nSOC )
    # for i in 1:size(data.Energy, 1)
    #     data.Energy[i,:] .*= data.k_w[i]
    # end
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    

    p = findfirst(x -> x > 0, data_energies)

    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 2000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    dos      = zeros(Float64, nbins)

    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            dos[j]    += gnorm * pdf(gauss, Ei - ebins[j])
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, dos, ebins, lw=2, label = false, lc = :black)
    plot!(p1,
        ylims=(-1.2,0.2),
        xlims=(0, maximum(dos) / 6),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid = false,
        xgrid = false, xticks = false,
        xlabel = L"g(E)"*" [states/eV]",
    )
    dos[1] = 0.0
    plot!(p1, dos[1:tmp], ebins[1:tmp], 
        fill = (0.2, -0.00000001, :gray),
        lw= LINE_WIDTH[1],
        lc=:black,
        label=false
        )
    return p1
end

function PDOS_PROCAR_nSOC()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.58199864
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_PDOS_nSOC( FILE_PROCAR_nSOC )
    for i in 1:size(data.Projected, 1)
        data.Projected[i,:] .*= data.k_w[i]
    end
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    data_projected  = vcat(data.Projected...)[p]

    
    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 2000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    pdos      = zeros(Float64, nbins)

    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            gval       = gnorm * pdf(gauss, Ei - ebins[j])
            pdos[j]   += gval * data_projected[idx]
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, pdos, ebins, lw=2, label = false)
    plot!(p1,
        ylims=(-1.2,0.2),
        xlims=(0, maximum(pdos) / 2 ),
        xtickfontsize  = 12,
        legendfontsize = 12,
        ygrid = false, yticks = false,
        xgrid = false, xticks = false,
        xlabel = "states/eV",
        title        = "PDOS"
    )
    pdos[1] = 0.0
    plot!(p1, pdos[1:tmp], ebins[1:tmp], 
        fill = (0.2, 0, :steelblue),
        lw= LINE_WIDTH[1],
        lc=1,
        label=false
        )
    return p1
end

####################################################################
######################## nSOC  + H #################################
####################################################################
" This function computes the DOS using the provided data obtained directly from the PROCAR file."
function DOS_PROCAR_nSOC_H()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.53063566
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_DOS_nSOC( FILE_PROCAR_nSOC_H )
    # for i in 1:size(data.Energy, 1)
    #     data.Energy[i,:] .*= data.k_w[i]
    # end
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    

    p = findfirst(x -> x > 0, data_energies)

    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 2000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    dos      = zeros(Float64, nbins)

    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            dos[j]    += gnorm * pdf(gauss, Ei - ebins[j])
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, dos, ebins, lw=2, label = false, lc = :black)
    plot!(p1,
        ylims=(-1.2,0.2),
        xlims=(0, maximum(dos) / 6),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid = false,
        xgrid = false, xticks = false,
        xlabel = L"g(E)"*" [states/eV]",
    )
    dos[1] = 0.0
    plot!(p1, dos[1:tmp], ebins[1:tmp], 
        fill = (0.2, -0.0000001, :gray),
        lw= LINE_WIDTH[1],
        lc=:black,
        label=false
        )
    return p1
end

function PDOS_PROCAR_nSOC_H()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.53063566
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data    = load_PROCAR_PDOS_nSOC_H( FILE_PROCAR_nSOC_H )
    for i in 1:size(data.Projected, 1)
        data.Projected[i,:] .*= data.k_w[i]
        data.ProjectedH[i,:] .*= data.k_w[i]
    end
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    data_projected  = vcat(data.Projected...)[p]
    data_projectedH = vcat(data.ProjectedH...)[p]

    
    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 2000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    pdos     = zeros(Float64, nbins)
    pdosH    = zeros(Float64, nbins)

    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            gval       = gnorm * pdf(gauss, Ei - ebins[j])
            pdos[j]   += gval * data_projected[idx]
            pdosH[j]  += gval * data_projectedH[idx]
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, ebins, pdos,  lw=2, lc = 1, label = false)
    plot!(p1, ebins, 10 .* pdosH, lw=2, lc = :red, label = false)
    plot!(p1,
        xlims=(-6.5,1),
        ylims=(-50, maximum(pdos) / 6 ),
        ytickfontsize  = 12,
        xtickfontsize  = 12,
        legendfontsize = 12,
        yguidefontsize=16,
        xguidefontsize=16,
        # ygrid = false, yticks = false,
        # xgrid = false, xticks = false,
        ylabel = "states/eV",
        xlabel = L"E-E_F"*" [eV]",
        title        = "non-SOC + H"
    )
    pdos[1] = 0.0
    plot!(p1, ebins[1:tmp], pdos[1:tmp],
        fill = (0.2, 0, :steelblue),
        lw= LINE_WIDTH[1],
        lc=1,
        label=false
        )
    plot!(p1, ebins[1:tmp], 10 .* pdosH[1:tmp],
        fill = (0.2, 0, :red),
        lw= LINE_WIDTH[1],
        lc=:red,
        label=false
        )
    ebins, pCOHP = plot_IpCOHP_nSOC()
    p = twinx()
    plot!(p, ebins, pCOHP, lw=2, lc=:green, label = false,
        xlims=(-6.5,1),ylabel = "-pCOHP", ytickfontsize  = 12,
        legendfontsize = 12, yguidefontsize=16, y_guidefontcolor=:green,
        y_foreground_color_axis=:green, y_foreground_color_text=:green,
        y_foreground_color_border=:green)

    return p1
end

###############################################################
######################### SOC #################################
###############################################################


" This function computes the DOS using the provided data obtained directly from the PROCAR file."
function load_PROCAR_DOS_SOC(path::String)
    num_kpoints = []; num_bnds = []; num_ions = []; Energies = []; k_weight = [];
    open(path, "r") do io
        readline(io)                            # comment / timestamp line
        tmp = split(readline(io))
        num_kpoints = parse(Int,tmp[4])
        num_bnds    = parse(Int,tmp[8])
        num_ions    = parse(Int,tmp[12])
        Energies    = zeros(num_kpoints, num_bnds)
        k_weight    = zeros(num_kpoints)
        readline(io)
        for i in 1:num_kpoints
            k_weight[i] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = split(readline(io))
                Energies[i, j] = parse(Float64, tmp[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                for i in 1:2+(1+num_ions)*4-1
                    readline(io)  # Skip the rest of the lines for this k-point
                end
            end
            readline(io)
        end
    end
    return PROCARDataDOS(num_kpoints, num_bnds, num_ions, Energies, k_weight)
end

function load_PROCAR_PDOS_SOC(path::String)
    num_kpoints = []; num_bnds = []; num_ions = [];
    k_weight = []; PDOS = []; Energies = [];
    open(path, "r") do io
        readline(io)                            # comment / timestamp line
        tmp = split(readline(io))
        num_kpoints = parse(Int,tmp[4])
        num_bnds    = parse(Int,tmp[8])
        num_ions    = parse(Int,tmp[12])
        Energies    = zeros(num_kpoints, num_bnds)
        PDOS        = zeros(num_kpoints, num_bnds)
        k_weight    = zeros(num_kpoints)
        readline(io)
        for i in 1:num_kpoints
            k_weight[i] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = readline(io)
                Energies[i, j] = parse(Float64, split(tmp)[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                tmp2 = 0
                for _ in 1:(1+num_ions)*4
                    tmp = readline(io)
                    # PDOS for the surface atoms
                    if split(tmp)[1] == "1" || split(tmp)[1] == "2" || split(tmp)[1] == "29" || split(tmp)[1] == "30"
                        if div(tmp2, 4) == 0  # s-orbital
                            PDOS[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        end
                        tmp2 += 1 
                    end
                end
                readline(io)  # Skip the rest of the lines for this k-point
            end
            readline(io)
        end
    end
    return PROCARDataPDOS(num_kpoints, num_bnds, num_ions, Energies, PDOS, k_weight)
end

function load_PROCAR_MAG_SOC(path::String)
    num_kpoints = []; num_bnds = []; num_ions = []; Energies = [];
    k_weight = []; PDOS = []; M_x = []; M_y = []; M_z = [];
    open(path, "r") do io
        readline(io)                            # comment / timestamp line
        tmp = split(readline(io))
        num_kpoints = parse(Int,tmp[4])
        num_bnds    = parse(Int,tmp[8])
        num_ions    = parse(Int,tmp[12])
        Energies    = zeros(num_kpoints, num_bnds)
        PDOS        = zeros(num_kpoints, num_bnds)
        M_x         = zeros(num_kpoints, num_bnds)
        M_y         = zeros(num_kpoints, num_bnds)
        M_z         = zeros(num_kpoints, num_bnds)
        k_weight    = zeros(num_kpoints)
        readline(io)
        for i in 1:num_kpoints
            k_weight[i] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = readline(io)
                Energies[i, j] = parse(Float64, split(tmp)[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                tmp2 = 0
                for _ in 1:(1+num_ions)*4
                    tmp = readline(io)
                    # PDOS for the surface atoms
                    if split(tmp)[1] == "1" || split(tmp)[1] == "2" || split(tmp)[1] == "29" || split(tmp)[1] == "30"
                        if div(tmp2, 4) == 0  # s-orbital
                            PDOS[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        elseif div(tmp2, 4) == 1  # p_x-orbital
                            M_x[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        elseif div(tmp2, 4) == 2  # p_y-orbital
                            M_y[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        elseif div(tmp2, 4) == 3  # p_z-orbital
                            M_z[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        end
                        tmp2 += 1 
                    end
                end
                readline(io)  # Skip the rest of the lines for this k-point
            end
            readline(io)
        end
    end
    return PROCARData_SOC_MAG(num_kpoints, num_bnds, num_ions, Energies, PDOS, M_x, M_y, M_z, k_weight)
end

function load_PROCAR_PDOS_SOC_H(path::String)
    num_kpoints = []; num_bnds = []; num_ions = [];
    Energies = []; k_weight = []; PDOS = []; PDOSH = []
    open(path, "r") do io
        readline(io)                            # comment / timestamp line
        tmp = split(readline(io))
        num_kpoints = parse(Int,tmp[4])
        num_bnds    = parse(Int,tmp[8])
        num_ions    = parse(Int,tmp[12])
        Energies    = zeros(num_kpoints, num_bnds)
        PDOS        = zeros(num_kpoints, num_bnds)
        PDOSH       = zeros(num_kpoints, num_bnds) 
        k_weight    = zeros(num_kpoints)
        readline(io)
        for i in 1:num_kpoints
            k_weight[i] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = readline(io)
                Energies[i, j] = parse(Float64, split(tmp)[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                tmp2  = 0
                tmp2H = 0
                for _ in 1:(1+num_ions)*4
                    tmp = readline(io)
                    # PDOS for the Hydrogen atom
                    if split(tmp)[1] == "1"
                        if tmp2H == 0  # s-orbital
                            PDOSH[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        end
                        tmp2H += 1
                    end 
                    # PDOS for the surface atoms
                    if split(tmp)[1] == "2" || split(tmp)[1] == "3" || split(tmp)[1] == "30" || split(tmp)[1] == "31"
                        if div(tmp2, 4) == 0  # Not the Magnetic moment
                            PDOS[i, j] += parse(Float64, split(tmp)[end])  # Read the PDOS value
                        end
                        tmp2 += 1 
                    end
                end
                readline(io)  # Skip the rest of the lines for this k-point
            end
            readline(io)
        end
    end
    return PROCARDataPDOSH(num_kpoints, num_bnds, num_ions, Energies, PDOS, PDOSH, k_weight)
end

function load_PROCAR_MAG_SOC_H(path::String)
    num_kpoints = []; num_bnds = []; num_ions = [];
    Energies = []; k_weight = []; PDOS = []; PDOSH = []
    M_x = []; M_y = []; M_z = [];
    open(path, "r") do io
        readline(io)                            # comment / timestamp line
        tmp = split(readline(io))
        num_kpoints = parse(Int,tmp[4])
        num_bnds    = parse(Int,tmp[8])
        num_ions    = parse(Int,tmp[12])
        Energies    = zeros(num_kpoints, num_bnds)
        PDOS        = zeros(num_kpoints, num_bnds)
        PDOSH       = zeros(num_kpoints, num_bnds) 
        M_x         = zeros(num_kpoints, num_bnds)
        M_y         = zeros(num_kpoints, num_bnds)
        M_z         = zeros(num_kpoints, num_bnds)
        k_weight    = zeros(num_kpoints)
        readline(io)
        for i in 1:num_kpoints
            k_weight[i] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = readline(io)
                Energies[i, j] = parse(Float64, split(tmp)[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                tmp2  = 0
                tmp2H = 0
                for _ in 1:(1+num_ions)*4
                    tmp = readline(io)
                    # PDOS for the Hydrogen atom
                    if split(tmp)[1] == "1"
                        if tmp2H == 0  # s-orbital
                            PDOSH[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        elseif tmp2H == 1  # σ_x-magnetization
                            M_x[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        elseif tmp2H == 2  # σ_y-magnetization
                            M_y[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        elseif tmp2H == 3  # σ_z-magnetization
                            M_z[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        end
                        tmp2H += 1
                    end 
                    # PDOS for the surface atoms
                    if split(tmp)[1] == "2" || split(tmp)[1] == "3" || split(tmp)[1] == "30" || split(tmp)[1] == "1"
                        if div(tmp2, 4) == 0  # s-orbital
                            PDOS[i, j] += parse(Float64, split(tmp)[end])  # Read the PDOS value
                        # elseif div(tmp2, 4) == 1  # σ_x-magnetization
                        #     M_x[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        # elseif div(tmp2, 4) == 2  # σ_y-magnetization
                        #     M_y[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        # elseif div(tmp2, 4) == 3  # σ_z-magnetization
                        #     M_z[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        end
                        tmp2 += 1 
                    end
                end
                readline(io)  # Skip the rest of the lines for this k-point
            end
            readline(io)
        end
    end
    return PROCARData_SOC_PDOS_H(num_kpoints, num_bnds, num_ions, Energies, PDOS, PDOSH, M_x, M_y, M_z, k_weight)
end

" This function computes the DOS using the provided data obtained directly from the PROCAR file."
function DOS_PROCAR_SOC()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.58485570
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_DOS_SOC( FILE_PROCAR_SOC )
    # for i in 1:size(data.Energy, 1)
    #     data.Energy[i,:] .*= data.k_w[i]
    # end
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    

    p = findfirst(x -> x > 0, data_energies)

    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 5000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    dos      = zeros(Float64, nbins)

    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            dos[j]    += gnorm * pdf(gauss, Ei - ebins[j])
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot(size = (200, 500))
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, dos, ebins, lw=2, label = false, lc = :black)
    plot!(p1,
        ylims=(-6.0,9.2),
        xlims=(0, maximum(dos)/4 ),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid = false,
        xgrid = false, xticks = false,
        xlabel = L"g(E)"*" [states/eV]",
    )
    dos[1] = 0.0
    plot!(p1, dos[1:tmp], ebins[1:tmp], 
        fill = (0.2, -0.000001, :gray),
        lw= LINE_WIDTH[1],
        lc=:black,
        label=false
        )
    return p1
end

function PDOS_PROCAR_SOC()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.58485570
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_PDOS_SOC( FILE_PROCAR_SOC )
    for i in 1:size(data.Projected, 1)
        data.Projected[i,:] .*= data.k_w[i]
    end
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    data_projected  = vcat(data.Projected...)[p]

    
    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 5000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    pdos      = zeros(Float64, nbins)

    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            gval       = gnorm * pdf(gauss, Ei - ebins[j])
            pdos[j]   += gval * data_projected[idx]
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot(size = (200, 500))
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, pdos, ebins, lw=2, label = false)
    plot!(p1,
        ylims=(-6.0,9.2),
        xlims=(0, maximum(pdos)/3 ),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid = false,
        xgrid = false, xticks = false,
        xlabel = L"g(E)"*" [states/eV]",
    )
    pdos[1] = 0.0
    plot!(p1, pdos[1:tmp], ebins[1:tmp], 
        fill = (0.2, 0, :steelblue),
        lw= LINE_WIDTH[1],
        lc=1,
        label=false
        )
    # plot!(p1, [0, pdos[tmp], pdos[tmp], 0.0], [-1.2, -1.2, ebins[tmp], ebins[tmp] ],
        # fillrange = 0.0,
        # fillalpha = 0.2,
        # fillcolor = :steelblue,
        # lw= 0,
        # lc=2,
        # label=false
        # )
    return p1
end

###############################################################
####################### SOC + H ###############################
###############################################################
" This function computes the DOS using the provided data obtained directly from the PROCAR file."
function DOS_PROCAR_SOC_H()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.53348380
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_DOS_SOC( FILE_PROCAR_SOC_H )
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    

    p = findfirst(x -> x > 0, data_energies)

    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 2000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    dos      = zeros(Float64, nbins)

    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            dos[j]    += gnorm * pdf(gauss, Ei - ebins[j])
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, dos, ebins, lw=2, label = false, lc = :black)
    plot!(p1,
        ylims=(-1.2,0.2),
        xlims=(0, maximum(dos) / 5),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid = false,
        xgrid = false, xticks = false,
        xlabel = L"g(E)"*" [states/eV]",
    )
    dos[1] = 0.0
    plot!(p1, dos[1:tmp], ebins[1:tmp], 
        fill = (0.2, -0.000000001, :gray),
        lw= LINE_WIDTH[1],
        lc=:black,
        label=fa1se
        )
    return p1
end

function PDOS_PROCAR_SOC_H()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.53063566
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_PDOS_SOC_H( FILE_PROCAR_SOC_H )
    for i in 1:size(data.Projected, 1)
        data.Projected[i,:]  .*= data.k_w[i]
        data.ProjectedH[i,:] .*= data.k_w[i]
    end
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    data_projected  = vcat(data.Projected...)[p]
    data_projectedH = vcat(data.ProjectedH...)[p]

    
    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 5000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    pdos     = zeros(Float64, nbins)
    pdosH    = zeros(Float64, nbins)
    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            gval       = gnorm * pdf(gauss, Ei - ebins[j])
            pdos[j]   += gval * data_projected[idx]
            pdosH[j]  += gval * data_projectedH[idx]
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot(size = (200, 500))
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, pdos, ebins, lw=2, label = false)
    plot!(p1, pdosH, ebins, lw=2, lc = :red, label = false)
    plot!(p1,
        ylims=(-7.0,5.0),
        xlims=(0, maximum(pdos)/4.5 ),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid = false,
        xgrid = false, xticks = false,
        xlabel = L"g(E)"*" [states/eV]",
    )
    pdos[1] = 0.0
    plot!(p1, pdos[1:tmp], ebins[1:tmp], 
        fill = (0.2, 0, :steelblue),
        lw= LINE_WIDTH[1],
        lc=1,
        label=false
        )
    plot!(p1, pdosH[1:tmp], ebins[1:tmp], 
        fill = (0.2, 0, :red),
        lw= LINE_WIDTH[1],
        lc=:red,
        label=false
        )
    return p1
end

function PDOS_PROCAR_SOC_H_pCOHP()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.53348380
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_PDOS_SOC_H( FILE_PROCAR_SOC_H )
    for i in 1:size(data.Projected, 1)
        data.Projected[i,:]  .*= data.k_w[i]
        data.ProjectedH[i,:] .*= data.k_w[i]
    end
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    data_projected  = vcat(data.Projected...)[p]
    data_projectedH = vcat(data.ProjectedH...)[p]

    
    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 2000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    pdos     = zeros(Float64, nbins)
    pdosH    = zeros(Float64, nbins)
    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            gval       = gnorm * pdf(gauss, Ei - ebins[j])
            pdos[j]   += gval * data_projected[idx]
            pdosH[j]  += gval * data_projectedH[idx]
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, ebins, pdos, lw=2, label = false)
    plot!(p1, ebins, 10 .* pdosH, lw=2, lc = :red, label = false)
    plot!(p1,
        framestyle=:box, 
        xlims=(-6.5,1),
        ylims=(-8, maximum(pdos) / 5.5 ),
        ytickfontsize  = 12,
        xtickfontsize  = 12,
        legendfontsize = 12,
        yguidefontsize=16,
        xguidefontsize=16,
        # ygrid = false, yticks = false,
        # xgrid = false, xticks = false,
        ylabel = "states/eV",
        xlabel = L"E-E_F"*" [eV]",
        title        = "SOC + H"
    )
    pdos[1] = 0.0
    plot!(p1, ebins[1:tmp], pdos[1:tmp],
        fill = (0.2, 0, :steelblue),
        lw= LINE_WIDTH[1],
        lc=1,
        label=false
        )
    plot!(p1, ebins[1:tmp], 10 .* pdosH[1:tmp], 
        fill = (0.2, 0, :red),
        lw= LINE_WIDTH[1],
        lc=:red,
        label=false
        )
    ebins, pCOHP = plot_IpCOHP_SOC()
    p = twinx()
    plot!(p, ebins, pCOHP, lw=2, lc=:green, label = false,
        xlims=(-6.5,1),ylabel = "-pCOHP", ytickfontsize  = 12,
        legendfontsize = 12, yguidefontsize=16, y_guidefontcolor=:green,
        y_foreground_color_axis=:green, y_foreground_color_text=:green,
        y_foreground_color_border=:green)
    return p1
end


" This function computes the DOS using the provided data obtained from VASPKIT."

function DOS_VASPKIT_SOC()
    data_tmp = readdlm(joinpath(PATH_FILES_H, "PDOS_A1_SOC.dat"))
    e_fermi     = 5.533484
    data_energy = data_tmp[2:end,1]
    data_s      = data_tmp[2:end,2]
    data_py     = data_tmp[2:end,3]
    data_pz     = data_tmp[2:end,4]
    data_px     = data_tmp[2:end,5]
    data_tot    = data_tmp[2:end,11]

    p = findfirst(x -> x > e_fermi, data_energy)
    println(sum(data_s[1:p-1]) / sum(data_tot[1:p-1]))

    ##############################################################################
    # 2. DOS-style smearing
    ##############################################################################
    # --- user parameters ---
    σ        = 0.09  # smearing width
    Emin     = minimum(data_energy)
    Emax     = maximum(data_energy)
    nbins    = 2000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    pdos_s   = zeros(Float64, nbins)
    pdos_py  = zeros(Float64, nbins)
    pdos_pz  = zeros(Float64, nbins)
    pdos_px  = zeros(Float64, nbins)
    for idx in 1:length(data_energy)
        Ei   = data_energy[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            gval       = gnorm * pdf(gauss, Ei - ebins[j])
            pdos_s[j]  += gval * data_s[idx]
            pdos_py[j] += gval * data_py[idx]
            pdos_pz[j] += gval * data_pz[idx]
            pdos_px[j] += gval * data_px[idx]
        end
    end
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:green, :blue, :black]
    plot!(p1, ebins, pdos_s,  label=L"H_{s}", lc = 1, ls = :solid)
    plot!(p1, ebins, pdos_px, label=L"H_{px}", lc = 2, ls = :dash)
    plot!(p1, ebins, pdos_py, label=L"H_{py}", lc = 3, ls = :dash)
    plot!(p1, ebins, pdos_pz, label=L"H_{pz}", lc = 4, ls = :dash)
    plot!(p1,
        xlims=(-10,10),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12
    )
    return p1
end

###############################################################
################# SOC + H Supercell 2x1x7 #####################
###############################################################

function load_PROCAR_PDOS_SOC_H_S(path::String)
    num_kpoints = []; num_bnds = []; num_ions = [];
    Energies = []; k_weight = []; PDOS = []; PDOSH = []
    open(path, "r") do io
        readline(io)                            # comment / timestamp line
        tmp = split(readline(io))
        num_kpoints = parse(Int,tmp[4])
        num_bnds    = parse(Int,tmp[8])
        num_ions    = parse(Int,tmp[12])
        Energies    = zeros(num_kpoints, num_bnds)
        PDOS        = zeros(num_kpoints, num_bnds)
        PDOSH       = zeros(num_kpoints, num_bnds) 
        k_weight    = zeros(num_kpoints)
        readline(io)
        for i in 1:num_kpoints
            k_weight[i] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = readline(io)
                Energies[i, j] = parse(Float64, split(tmp)[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                tmp2  = 0
                tmp2H = 0
                for _ in 1:(1+num_ions)*4
                    tmp = readline(io)
                    # PDOS for the Hydrogen atom
                    if split(tmp)[1] == "1"
                        if tmp2H == 0  # s-orbital
                            PDOSH[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        end
                        tmp2H += 1
                    end 
                    # PDOS for the surface atoms
                    if split(tmp)[1] == "2" || split(tmp)[1] == "3" || split(tmp)[1] == "4" || split(tmp)[1] == "5" || split(tmp)[1] == "58" || split(tmp)[1] == "59" || split(tmp)[1] == "60" || split(tmp)[1] == "61"
                        if div(tmp2, 8) == 0  # s-orbital
                            PDOS[i, j] += parse(Float64, split(tmp)[end])  # Read the PDOS value
                        end
                        tmp2 += 1 
                    end
                end
                readline(io)  # Skip the rest of the lines for this k-point
            end
            readline(io)
        end
    end
    return PROCARDataPDOSH(num_kpoints, num_bnds, num_ions, Energies, PDOS, PDOSH, k_weight)
end

function load_PROCAR_MAG_SOC_H_S(path::String)
    num_kpoints = []; num_bnds = []; num_ions = [];
    Energies = []; k_weight = []; PDOS = []; PDOSH = []
    M_x = []; M_y = []; M_z = [];
    open(path, "r") do io
        readline(io)                            # comment / timestamp line
        tmp = split(readline(io))
        num_kpoints = parse(Int,tmp[4])
        num_bnds    = parse(Int,tmp[8])
        num_ions    = parse(Int,tmp[12])
        Energies    = zeros(num_kpoints, num_bnds)
        PDOS        = zeros(num_kpoints, num_bnds)
        PDOSH       = zeros(num_kpoints, num_bnds) 
        M_x         = zeros(num_kpoints, num_bnds)
        M_y         = zeros(num_kpoints, num_bnds)
        M_z         = zeros(num_kpoints, num_bnds)
        k_weight    = zeros(num_kpoints)
        readline(io)
        for i in 1:num_kpoints
            k_weight[i] = parse(Float64, split(readline(io))[end])  # Skip the k-point line
            readline(io)  # Skip the comment line
            for j in 1:num_bnds
                tmp = readline(io)
                Energies[i, j] = parse(Float64, split(tmp)[5])  # Read the energy value
                readline(io)  # Skip the line with the band index
                readline(io)  # Skip the line with the band index
                tmp2  = 0
                tmp2H = 0
                for _ in 1:(1+num_ions)*4
                    tmp = readline(io)
                    # PDOS for the Hydrogen atom
                    if split(tmp)[1] == "1"
                        if tmp2H == 0  # s-orbital
                            PDOSH[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        elseif tmp2H == 1  # σ_x-magnetization
                            M_x[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        elseif tmp2H == 2  # σ_y-magnetization
                            M_y[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        elseif tmp2H == 3  # σ_z-magnetization
                            M_z[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        end
                        tmp2H += 1
                    end 
                    # PDOS for the surface atoms
                    if split(tmp)[1] == "2" || split(tmp)[1] == "3" || split(tmp)[1] == "30" || split(tmp)[1] == "1"
                        if div(tmp2, 4) == 0  # s-orbital
                            PDOS[i, j] += parse(Float64, split(tmp)[end])  # Read the PDOS value
                        # elseif div(tmp2, 4) == 1  # σ_x-magnetization
                        #     M_x[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        # elseif div(tmp2, 4) == 2  # σ_y-magnetization
                        #     M_y[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        # elseif div(tmp2, 4) == 3  # σ_z-magnetization
                        #     M_z[i, j] += parse(Float64, split(tmp)[end])  # Read the energy value
                        end
                        tmp2 += 1 
                    end
                end
                readline(io)  # Skip the rest of the lines for this k-point
            end
            readline(io)
        end
    end
    return PROCARData_SOC_PDOS_H(num_kpoints, num_bnds, num_ions, Energies, PDOS, PDOSH, M_x, M_y, M_z, k_weight)
end

function DOS_PROCAR_SOC_H_S()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.74715706
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_DOS_SOC( FILE_PROCAR_SOC_H_S )
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    

    p = findfirst(x -> x > 0, data_energies)

    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 2000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    dos      = zeros(Float64, nbins)

    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            dos[j]    += gnorm * pdf(gauss, Ei - ebins[j])
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, dos, ebins, lw=2, label = false, lc = :black)
    plot!(p1,
        ylims=(-1.2,0.2),
        xlims=(0, maximum(dos) / 5),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 18,
        legendfontsize = 12,
        ygrid = false,
        xgrid = false, xticks = false,
        xlabel = L"g(E)"*" [states/eV]",
    )
    dos[1] = 0.0
    plot!(p1, dos[1:tmp], ebins[1:tmp], 
        fill = (0.2, -0.0000001, :gray),
        lw= LINE_WIDTH[1],
        lc= :black,
        label=false
        )
    return p1
end

function PDOS_PROCAR_SOC_H_S()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.74715706
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_PDOS_SOC_H_S( FILE_PROCAR_SOC_H_S )
    for i in 1:size(data.Projected, 1)
        data.Projected[i,:]  .*= data.k_w[i]
        data.ProjectedH[i,:] .*= data.k_w[i]
    end
    p               = sortperm( vcat( data.Energy... ) )
    data_energies   = vcat(data.Energy...)[p] .- e_fermi
    data_projected  = vcat(data.Projected...)[p]
    data_projectedH = vcat(data.ProjectedH...)[p]

    
    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    σ        = 0.04  # smearing width
    Emin     = minimum(data_energies)
    Emax     = maximum(data_energies)
    nbins    = 2000  # number of bins
    # ----------------------
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    pdos     = zeros(Float64, nbins)
    pdosH    = zeros(Float64, nbins)
    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            gval       = gnorm * pdf(gauss, Ei - ebins[j])
            pdos[j]   += gval * data_projected[idx]
            pdosH[j]  += gval * data_projectedH[idx]
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, pdos, ebins, lw=2, label = false)
    plot!(p1, pdosH, ebins, lw=2, lc = :red, label = false)
    plot!(p1,
        ylims=(-1.2,0.2),
        xlims=(0, maximum(pdos) / 5.0 ),
        xtickfontsize  = 12,
        legendfontsize = 12,
        ygrid = false, yticks = false,
        xgrid = false, xticks = false,
        xlabel = "states/eV",
        title        = "PDOS"
    )
    pdos[1] = 0.0
    plot!(p1, pdos[1:tmp], ebins[1:tmp], 
        fill = (0.2, 0, :steelblue),
        lw= LINE_WIDTH[1],
        lc=1,
        label=false
        )
    plot!(p1, pdosH[1:tmp], ebins[1:tmp], 
        fill = (0.2, 0, :red),
        lw= LINE_WIDTH[1],
        lc=:red,
        label=false
        )
    return p1
end

function tst()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.53348380
    T       = 0.0

    ####################################################
    # 1. Load the PROCAR file
    ####################################################
    data            = load_PROCAR_PDOS_SOC_H( FILE_PROCAR_SOC_H )
    for i in 1:size(data.Projected, 1)
        data.Projected[i,:]  .*= data.k_w[i]
        data.ProjectedH[i,:] .*= data.k_w[i]
    end

    return sum( FD_dist.(data.Energy,0.0,e_fermi) .* data.ProjectedH )
end