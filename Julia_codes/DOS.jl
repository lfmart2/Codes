using LaTeXStrings
using Plots
using Distributions
using DelimitedFiles
using HDF5

const FILE_DATA_nSOC      = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/Collinear/00_Data"
const FILE_DATA_nSOC_H    = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/Collinear_H/00_Data"
const FILE_DATA_SOC       = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/nCollinear/00_Data"
const FILE_DATA_SOC_H     = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/nCollinear_H/00_Data"
const FILE_DATA_SOC_H_S   = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell_217/nCollinear_H/00_Data"

####################################################################
######################### nSOC #####################################
####################################################################
" This function computes the DOS using the provided data obtained directly from the PROCAR file."
function DOS_nSOC()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.645736
    T       = 0.0
    k       = 1:253
    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_nSOC_H,"FS_H_up_201x201.h5"), "r")
    Energies = read(fid, "eigenvalues")[:,:,k] 
    Energies .= Energies .- e_fermi

    close(fid)
    # for i in 1:size(data.Energy, 1)
    #     data.Energy[i,:] .*= data.k_w[i]
    # end
    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi
    

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
        xlims=(0, maximum(dos)./ 1.5),
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
    e_fermi = 5.645736
    T       = 0.0
    k       = 1:253
    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_nSOC_H,"FS_H_up_201x201.h5"), "r")
    Energies = read(fid, "eigenvalues")[:,:,k] 
    Energies .= Energies .- e_fermi

    # Initialize arrays to store reconstructed data
    proj0_H = read(fid, "proj_H_up")[:,:,k]

    proj   = read(fid, "proj1")[:,:,k] .+ read(fid, "proj2")[:,:,k] .+ read(fid, "proj3")[:,:,k] .+ read(fid, "proj4")[:,:,k]

    close(fid)

    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi
    proj0_H  = vcat(proj0_H...)[p]
    proj  = vcat(proj...)[p]

    
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
            pdos[j]   += gval * proj[idx]
            pdosH[j]  += gval * proj0_H[idx]

        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, pdos, ebins, lw=2, label = false, lc = 1)
    plot!(p1, pdosH, ebins, lw=2, label = false, lc = :red)
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
    pdosH[1] = 0.0
    plot!(p1, pdosH[1:tmp], ebins[1:tmp], 
        fill = (0.2, 0, :red),
        lw= LINE_WIDTH[1],
        lc=:red,
        label=false
        )
    return p1
end

####################################################################
######################## nSOC  + H #################################
####################################################################
" This function computes the DOS using the provided data obtained directly from the PROCAR file."
function DOS_nSOC(Emin::Float64, Emax::Float64, nbins::Int64, σ::Float64; e_fermi::Float64 = 5.53063566)
    ####################################################
    # 0. Definitions
    ####################################################
    T       = 0.0
    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_nSOC_H,"DOS_H_up_201x201.h5"), "r")
    Energies_up = read(fid, "eigenvalues") .- e_fermi

    close(fid)

    fid = h5open(joinpath(FILE_DATA_nSOC_H,"DOS_H_dn_201x201.h5"), "r")
    Energies_dn = read(fid, "eigenvalues") .- e_fermi

    close(fid)
    # for i in 1:size(data.Energy, 1)
    #     data.Energy[i,:] .*= data.k_w[i]
    # end
    p_up             = sortperm( vcat( Energies_up... ) )
    p_dn             = sortperm( vcat( Energies_dn... ) )
    data_energies_up = vcat(Energies_up...)[p_up]
    data_energies_dn = vcat(Energies_dn...)[p_dn]
    


    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- calculate the smearing ---
    ebins    = range(Emin, Emax, length=nbins)
    ΔE       = step(ebins)
    gnorm    = 1 / (σ * sqrt(2 * π))
    gauss    = Normal(0, σ)
    dos_up   = zeros(Float64, nbins)
    dos_dn   = zeros(Float64, nbins)

    for idx in 1:length(data_energies_up)
        Ei_up = data_energies_up[idx]
        Ei_dn = data_energies_dn[idx]
        # find the bin index for the current spin-up energy
        jmin_up = max(1, floor(Int, (Ei_up - 4 * σ - Emin) / ΔE) + 1)
        jmax_up = min(nbins, ceil(Int, (Ei_up + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin_up:jmax_up
            dos_up[j] += gnorm * pdf(gauss, Ei_up - ebins[j])
        end
        # find the bin index for the current spin-up energy
        jmin_dn = max(1, floor(Int, (Ei_dn - 4 * σ - Emin) / ΔE) + 1)
        jmax_dn = min(nbins, ceil(Int, (Ei_dn + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin_dn:jmax_dn
            dos_dn[j] += gnorm * pdf(gauss, Ei_dn - ebins[j])
        end
    end
    
    return dos_up, dos_dn
end

function DOS_PROCAR_nSOC_H()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.645736
    T       = 0.0
    k       = 1:253
    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_nSOC_H,"FS_H_up_201x201.h5"), "r")
    Energies = read(fid, "eigenvalues")[:,:,k] 
    Energies .= Energies .- e_fermi

    close(fid)
    # for i in 1:size(data.Energy, 1)
    #     data.Energy[i,:] .*= data.k_w[i]
    # end
    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi
    

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
    k       = 1:253
    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_nSOC_H,"FS_H_up_201x201.h5"), "r")
    Energies = read(fid, "eigenvalues")[:,:,k] 
    Energies .= Energies .- e_fermi

    # Initialize arrays to store reconstructed data
    proj0_H = read(fid, "proj_H_up")[:,:,k]

    proj   = read(fid, "proj1")[:,:,k] .+ read(fid, "proj2")[:,:,k] .+ read(fid, "proj3")[:,:,k] .+ read(fid, "proj4")[:,:,k]

    close(fid)

    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi
    proj0_H  = vcat(proj0_H...)[p]
    proj  = vcat(proj...)[p]

    
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
            pdos[j]   += gval * proj[idx]
            pdosH[j]  += gval * proj0_H[idx]

        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, pdos, ebins, lw=2, label = false, lc = 1)
    plot!(p1, pdosH, ebins, lw=2, label = false, lc = :red)
    plot!(p1,
        ylims=(-1.2,0.2),
        xlims=(0, maximum(pdos) / 1.5 ),
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
    pdosH[1] = 0.0
    plot!(p1, pdosH[1:tmp], ebins[1:tmp], 
        fill = (0.2, 0, :red),
        lw= LINE_WIDTH[1],
        lc=:red,
        label=false
        )
    return p1
end

###############################################################
######################### SOC #################################
###############################################################


" This function computes the DOS using the provided data obtained directly from the PROCAR file."
function DOS_SOC(Emin::Float64, Emax::Float64, nbins::Int64, σ::Float64; e_fermi::Float64 = 5.53063566)
    ####################################################
    # 0. Definitions
    ####################################################
    T       = 0.0

    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_SOC,"DOS_101x101.h5"), "r")
    Energies = read(fid, "eigenvalues") .- e_fermi

    close(fid)

    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi

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
        xlims=(0, maximum(dos)/1.5),
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

function DOS_PROCAR_SOC()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.62997313
    T       = 0.0
    k       = 1:504
    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_SOC,"DOS_101x101.h5"), "r")
    Energies = read(fid, "eigenvalues")[:,:,k] 
    Energies .= Energies .- e_fermi

    close(fid)

    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi

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
        xlims=(0, maximum(dos)/1.5),
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
    e_fermi = 5.62997313
    T       = 0.0
    k       = 1:504
    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_SOC,"DOS_101x101.h5"), "r")
    Energies = read(fid, "eigenvalues")[:,:,k] 
    Energies .= Energies .- e_fermi

    # Initialize arrays to store reconstructed data
    proj   = read(fid, "proj1")[:,:,k] .+ read(fid, "proj2")[:,:,k] .+ read(fid, "proj3")[:,:,k] .+ read(fid, "proj4")[:,:,k]

    close(fid)

    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi
    proj  = vcat(proj...)[p]

    
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

    for idx in 1:length(data_energies)
        Ei   = data_energies[idx]
        # find the bin index for the current energy
        jmin = max(1, floor(Int, (Ei - 4 * σ - Emin) / ΔE) + 1)
        jmax = min(nbins, ceil(Int, (Ei + 4 * σ - Emin) / ΔE) + 1)
        for j in jmin:jmax
            gval       = gnorm * pdf(gauss, Ei - ebins[j])
            pdos[j]   += gval * proj[idx]
        end
    end
    tmp = findmin(abs.(ebins))[2]


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, pdos, ebins, lw=2, label = false, lc = 1)
    plot!(p1,
        ylims=(-1.2,0.2),
        xlims=(0, maximum(pdos) / 1.2 ),
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

###############################################################
####################### SOC + H ###############################
###############################################################
" This function computes the DOS using the provided data obtained directly from the PROCAR file."
function DOS_SOC_H(Emin::Float64, Emax::Float64, nbins::Int64, σ::Float64; e_fermi::Float64 = 5.53348380)
    ####################################################
    # 0. Definitions
    ####################################################
    T       = 0.0

    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_SOC_H,"DOS_H_101x101.h5"), "r")
    Energies = read(fid, "eigenvalues") .- e_fermi

    close(fid)
    # for i in 1:size(data.Energy, 1)
    #     data.Energy[i,:] .*= data.k_w[i]
    # end
    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi
    
    ####################################################
    # 2. DOS-style smearing
    ####################################################
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

    return dos

end

function plot_DOS_SOC_H()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.53348380
    T       = 0.0
    k       = 1:506

    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_SOC_H,"DOS_H_101x101.h5"), "r")
    Energies = read(fid, "eigenvalues")[:,:,k] 
    Energies .= Energies .- e_fermi

    close(fid)
    # for i in 1:size(data.Energy, 1)
    #     data.Energy[i,:] .*= data.k_w[i]
    # end
    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi
    
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
        xlims=(0, maximum(dos) / 1.2),
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
        label = false
        )
    return p1
end

function PDOS_PROCAR_SOC_H(Emin::Float64, Emax::Float64, nbins::Int64, σ::Float64)
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.53348380
    T       = 0.0
    k       = 1:506
    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_SOC_H,"DOS_H_101x101.h5"), "r")
    Energies = read(fid, "eigenvalues")[:,:,k] 
    Energies .= Energies .- e_fermi

    # Initialize arrays to store reconstructed data
    proj0_H = read(fid, "proj_H_up")[:,:,k] .+ read(fid, "proj_H_dn")[:,:,k]

    proj   = read(fid, "proj1")[:,:,k] .+ read(fid, "proj2")[:,:,k] .+ read(fid, "proj3")[:,:,k] .+ read(fid, "proj4")[:,:,k]

    close(fid)
    
    nkx,nky,nbnds = size(Energies)

    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi
    proj0_H  = vcat(proj0_H...)[p]
    proj  = vcat(proj...)[p]

    
    ####################################################
    # 2. DOS-style smearing
    ####################################################
    # --- user parameters ---
    # σ        = 0.04  # smearing width
    # Emin     = minimum(data_energies)
    # Emax     = maximum(data_energies)
    # nbins    = 2000  # number of bins
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
            pdos[j]   += gval * proj[idx]
            pdosH[j]  += gval * proj0_H[idx]
        end
    end
    tmp = findmin(abs.(ebins))[2]

    pdos .= pdos ./ (nkx * nky)
    pdosH .= pdosH ./ (nkx * nky)


    ####################################################
    # 3. Plot the results
    ####################################################
    p1 = plot()
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    plot!(p1, pdos, ebins, lw=2, label = false, lc = 1)
    plot!(p1, pdosH, ebins, lw=2, label = false, lc = :red)
    plot!(p1,
        ylims=(-1.2,0.2),
        xlims=(0, maximum(pdos) / 1.2 ),
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
    pdosH[1] = 0.0
    plot!(p1, pdosH[1:tmp], ebins[1:tmp], 
        fill = (0.2, 0, :red),
        lw= LINE_WIDTH[1],
        lc=:red,
        label=false
        )
    return p1, pdos, pdosH, ebins
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

function DOS_PROCAR_SOC_H_S()
    ####################################################
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.74715706
    T       = 0.0

    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_SOC_H_S,"DOS_H_51x101.h5"), "r")
    Energies = read(fid, "eigenvalues")[:,:,k] 
    Energies .= Energies .- e_fermi

    close(fid)

    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi

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
        xlims=(0, maximum(dos) / 1.2),
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
    # 0. Definitions
    ####################################################
    # e_fermi = 5.533484 
    e_fermi = 5.53348380
    T       = 0.0
    k       = 500:1010
    ####################################################
    # 1. Load the .h5 file
    ####################################################
    fid = h5open(joinpath(FILE_DATA_SOC_H_S,"DOS_H_51x101.h5"), "r")
    Energies = read(fid, "eigenvalues")[:,:,k] 
    Energies .= Energies .- e_fermi

    # Initialize arrays to store reconstructed data
    proj0_H = read(fid, "proj_H_up")[:,:,k] .+ read(fid, "proj_H_dn")[:,:,k]

    proj   = read(fid, "proj1")[:,:,k] .+ read(fid, "proj2")[:,:,k] .+ read(fid, "proj3")[:,:,k] .+ read(fid, "proj4")[:,:,k]

    close(fid)

    p               = sortperm( vcat( Energies... ) )
    data_energies   = vcat(Energies...)[p] .- e_fermi
    proj0_H  = vcat(proj0_H...)[p]
    proj  = vcat(proj...)[p]

    
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
            pdos[j]   += gval * proj[idx]
            pdosH[j]  += gval * proj0_H[idx]
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
        xlims=(0, maximum(pdos) / 1.2 ),
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