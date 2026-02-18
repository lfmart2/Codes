using LaTeXStrings
using Plots, Contour
using Distributions
using DelimitedFiles
using HDF5

const FILE_DATA_nSOC      = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/Collinear/00_Data"
const FILE_DATA_nSOC_H    = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/Collinear_H/00_Data"
const FILE_DATA_SOC       = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/nCollinear/00_Data"
const FILE_DATA_SOC_H     = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell/nCollinear_bridge_gpu/00_Data"
const FILE_DATA_SOC_H_S   = raw"/media/r_floren/FernandoM/VASP/NbP/Supercell_217/nCollinear_H/00_Data"

""" This function computes the simple Fermi surface lines using the provided data obtained directly from the h5 files. Here, we have painted the 
    layer projection weights on top of the Fermi surface lines. The gray lines are zero energy modes in the bulk. """

####################################################################
######################### nSOC #####################################
####################################################################

function plot_FS_collinear(kstart::Int=50, kend::Int=252, ΔE_F::Float64=0.0)
    # ------------------------------------------------------------
    # 1.  Read the slice we need
    # ------------------------------------------------------------
    kslice = kstart:kend 
    E_F    = 5.58199864 + ΔE_F
    h5f    = h5open(joinpath(FILE_DATA_nSOC, "FS_up_201x201.h5"),"r")
    E_raw_up  = read(h5f,"eigenvalues_up")[:,:,kslice]
    Energies_up = real.(E_raw_up .- E_F)

    layer_projs_up = real.(sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4))
    close(h5f)

    h5f         = h5open(joinpath(FILE_DATA_nSOC, "FS_dn_201x201.h5"),"r")
    E_raw_dn    = read(h5f,"eigenvalues_dn")[:,:,kslice]
    Energies_dn = real.(E_raw_dn .- E_F)

    layer_projs_dn = real.(sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4))

    close(h5f)

    # sizes
    nkx, nky, nbands = size(Energies_up)

    # ------------------------------------------------------------
    # 2.  Mirror the first quadrant → full BZ helper
    # ------------------------------------------------------------
    mirror_q1(arr) = [reverse(arr)         reverse(arr; dims=1);
                      reverse(arr; dims=2) arr]

    E_full_up = cat((mirror_q1(@view Energies_up[:,:,b]) for b in 1:nbands)...; dims=3)
    E_full_dn = cat((mirror_q1(@view Energies_dn[:,:,b]) for b in 1:nbands)...; dims=3)
    W_full_up = cat((mirror_q1(@view layer_projs_up[:,:,b]) for b in 1:nbands)...; dims=3)
    W_full_dn = cat((mirror_q1(@view layer_projs_dn[:,:,b]) for b in 1:nbands)...; dims=3)

    # k‑grid
    kx = range(-π, π; length=size(E_full_up,1))
    ky = range(-π, π; length=size(E_full_up,2))
    Δk = step(kx)

    # ------------------------------------------------------------
    # 3.  Small helper to map a coordinate → nearest grid index
    # ------------------------------------------------------------
    idx(x, g) = clamp(Int(round((x - first(g)) / Δk)) + 1, 1, length(g))

    # ------------------------------------------------------------
    # 4.  Build the plot
    # ------------------------------------------------------------
    plt = plot(; xlabel=L"\frac{k_x}{a}", ylabel = L"\frac{k_y}{a}",
                #  xlims = (-4*π/5 + π, π), ylims = (0,1),
                 ylims = (-6*π/5, 6*π/5), xlims = (-6*π/5, 6*π/5),
                 xticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                 yticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                 framestyle=:box, aspect_ratio=:equal,
                 grid=false, size=(600,600))

    for b in 1:nbands
        for cl in levels(contours(kx,ky,E_full_up[:,:,b],[0.0]))
            for ln in lines(cl)
                xs, ys = coordinates(ln)
                vs     = [ W_full_up[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]

                plot!(plt, xs .- 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)

                plot!(plt, xs .- 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)


                # faint guide line
                plot!(plt, xs, ys; lc=:black, lw=1, alpha=0.15, label=false)
                # weight‑coloured overlay
                plot!(plt, xs, ys; lc=1, lw=2, linealpha=vs, label=false)
                

                plot!(plt, xs, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                

                plot!(plt, xs, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                

                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)

                plot!(plt, xs .+ 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                

                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
            end
        end
    end
    # Plot for Spin Down
    for b in 1:nbands
        for cl in levels(contours(kx,ky,E_full_dn[:,:,b],[0.0]))
            for ln in lines(cl)
                xs, ys = coordinates(ln)
                vs     = [ W_full_dn[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]

                plot!(plt, xs .- 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)

                plot!(plt, xs .- 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)


                # faint guide line
                plot!(plt, xs, ys; lc=:black, lw=1, alpha=0.15, label=false)
                # weight‑coloured overlay
                plot!(plt, xs, ys; lc=1, lw=2, linealpha=vs, label=false)

                plot!(plt, xs, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                
                plot!(plt, xs, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                

                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)

                plot!(plt, xs .+ 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
            end
        end
    end

    plot!(plt, [-pi; -pi], [pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [-pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi; pi], [-pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [0; π], [0; π], lc = :black, lw = 1, label = false, ls = :dashdot)
    # savefig(plt, "/media/r_floren/FernandoM/VASP/NbP/Supercell/Figures/Fermi_Surface/New_Calculations/Slab_H/TB/FS_surf.png")
    return plt
end


####################################################################
######################### nSOC + H #################################
####################################################################

function plot_FS_H_collinear(kstart::Int=50, kend::Int=253, ΔE_F::Float64=0.0)
    # ------------------------------------------------------------
    # 1.  Read the slice we need
    # ------------------------------------------------------------
    kslice = kstart:kend 
    E_F    = 5.53063566 + ΔE_F
    h5f    = h5open(joinpath(FILE_DATA_nSOC_H, "FS_H_up_201x201.h5"),"r")
    E_raw_up  = read(h5f,"eigenvalues")[:,:,kslice]
    Energies_up = real.(E_raw_up .- E_F)

    layer_projs_up = real.(sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4))
    Hlayer_proj_up = read(h5f,"proj_H_up")[:,:,kslice]
    close(h5f)

    h5f         = h5open(joinpath(FILE_DATA_nSOC_H, "FS_H_dn_201x201.h5"),"r")
    E_raw_dn    = read(h5f,"eigenvalues")[:,:,kslice]
    Energies_dn = real.(E_raw_dn .- E_F)

    layer_projs_dn = real.(sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4))
    Hlayer_proj_dn = read(h5f,"proj_H_dn")[:,:,kslice]

    close(h5f)

    # sizes
    nkx, nky, nbands = size(Energies_up)

    # ------------------------------------------------------------
    # 2.  Mirror the first quadrant → full BZ helper
    # ------------------------------------------------------------
    mirror_q1(arr) = [reverse(arr)         reverse(arr; dims=1);
                      reverse(arr; dims=2) arr]

    E_full_up = cat((mirror_q1(@view Energies_up[:,:,b]) for b in 1:nbands)...; dims=3)
    E_full_dn = cat((mirror_q1(@view Energies_dn[:,:,b]) for b in 1:nbands)...; dims=3)
    W_full_up = cat((mirror_q1(@view layer_projs_up[:,:,b]) for b in 1:nbands)...; dims=3)
    W_full_dn = cat((mirror_q1(@view layer_projs_dn[:,:,b]) for b in 1:nbands)...; dims=3)

    H_full_up = cat((mirror_q1(@view Hlayer_proj_up[:,:,b]) for b in 1:nbands)...; dims=3) .* 20
    H_full_dn = cat((mirror_q1(@view Hlayer_proj_dn[:,:,b]) for b in 1:nbands)...; dims=3) .* 20

    # k‑grid
    kx = range(-π, π; length=size(E_full_up,1))
    ky = range(-π, π; length=size(E_full_up,2))
    Δk = step(kx)

    # ------------------------------------------------------------
    # 3.  Small helper to map a coordinate → nearest grid index
    # ------------------------------------------------------------
    idx(x, g) = clamp(Int(round((x - first(g)) / Δk)) + 1, 1, length(g))

    # ------------------------------------------------------------
    # 4.  Build the plot
    # ------------------------------------------------------------
    plt = plot(; xlabel=L"\frac{k_x}{a}", ylabel = L"\frac{k_y}{a}",
                #  xlims = (-4*π/5 + π, π), ylims = (0,1),
                 ylims = (-6*π/5, 6*π/5), xlims = (-6*π/5, 6*π/5),
                 xticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                 yticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                 framestyle=:box, aspect_ratio=:equal,
                 grid=false, size=(600,600))

    for b in 1:nbands
        for cl in levels(contours(kx,ky,E_full_up[:,:,b],[0.0]))
            for ln in lines(cl)
                xs, ys = coordinates(ln)
                vs     = [ W_full_up[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]
                vs_H   = [ H_full_up[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]

                plot!(plt, xs .- 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)

                plot!(plt, xs .- 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- 2*π, ys; lc=2, lw=2, linealpha=vs_H, label=false)
                
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)


                # faint guide line
                plot!(plt, xs, ys; lc=:black, lw=1, alpha=0.15, label=false)
                # weight‑coloured overlay
                plot!(plt, xs, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys; lc=2, lw=2, linealpha=vs_H, label=false)
                

                plot!(plt, xs, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)
                

                plot!(plt, xs, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)
                

                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)

                plot!(plt, xs .+ 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=2, lw=2, linealpha=vs_H, label=false)
                

                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)
            end
        end
    end
    # Plot for Spin Down
    for b in 1:nbands
        for cl in levels(contours(kx,ky,E_full_dn[:,:,b],[0.0]))
            for ln in lines(cl)
                xs, ys = coordinates(ln)
                vs     = [ W_full_dn[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]
                vs_H   = [ H_full_dn[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]

                plot!(plt, xs .- 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)

                plot!(plt, xs .- 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- 2*π, ys; lc=2, lw=2, linealpha=vs_H, label=false)
                
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)


                # faint guide line
                plot!(plt, xs, ys; lc=:black, lw=1, alpha=0.15, label=false)
                # weight‑coloured overlay
                plot!(plt, xs, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys; lc=2, lw=2, linealpha=vs_H, label=false)

                plot!(plt, xs, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)
                
                plot!(plt, xs, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)
                

                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)

                plot!(plt, xs .+ 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=2, lw=2, linealpha=vs_H, label=false)
                
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)
            end
        end
    end

    plot!(plt, [-pi; -pi], [pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [-pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi; pi], [-pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [0; π], [0; π], lc = :black, lw = 1, label = false, ls = :dashdot)
    # savefig(plt, "/media/r_floren/FernandoM/VASP/NbP/Supercell/Figures/Fermi_Surface/New_Calculations/Slab_H/TB/FS_surf.png")
    return plt
end

####################################################################
############################## SOC #################################
####################################################################

" This function computes the DOS using the provided data obtained directly from the PROCAR file."

function plot_FS_ncollinear(kstart::Int=200, kend::Int=300, ΔE_F::Float64=0.0)
    # ------------------------------------------------------------
    # 1.  Read the slice we need
    # ------------------------------------------------------------
    kslice = kstart:kend
    h5f    = h5open(joinpath(FILE_DATA_SOC, "FS_201x201.h5"),"r")
    E_raw  = read(h5f,"eigenvalues")[:,:,kslice]
    E_F    = 5.62997313 - ΔE_F
    Energies = E_raw .- E_F

    layer_projs = sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4)
    # layer_projs = read(h5f, "proj1_dxy_up")[:,:,kslice] .+ read(h5f, "proj1_dxy_dn")[:,:,kslice]

    


    close(h5f)

    # sizes
    nkx, nky, nbands = size(Energies)

    # ------------------------------------------------------------
    # 2.  Mirror the first quadrant → full BZ helper
    # ------------------------------------------------------------
    mirror_q1(arr) = [reverse(arr)         reverse(arr; dims=1);
                      reverse(arr; dims=2) arr]

    E_full = cat((mirror_q1(@view Energies[:,:,b]) for b in 1:nbands)...; dims=3)
    W_full = cat((mirror_q1(@view layer_projs[:,:,b]) for b in 1:nbands)...; dims=3)

    # k‑grid
    kx = range(-π, π; length=size(E_full,1))
    ky = range(-π, π; length=size(E_full,2))
    Δk = step(kx)

    # ------------------------------------------------------------
    # 3.  Small helper to map a coordinate → nearest grid index
    # ------------------------------------------------------------
    idx(x, g) = clamp(Int(round((x - first(g)) / Δk)) + 1, 1, length(g))

    # ------------------------------------------------------------
    # 4.  Build the plot
    # ------------------------------------------------------------
    plt = plot(; xlabel=L"\frac{k_x}{a}", ylabel = L"\frac{k_y}{a}",
                 xlims = (-6*π/5,6*π/5), ylims = (-6*π/5,6*π/5),
                 xticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                 yticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                framestyle=:box, aspect_ratio=:equal,
                 grid=false, size=(600,600))

    for b in 1:nbands
        for cl in levels(contours(kx,ky,E_full[:,:,b],[-0.00]))
            for ln in lines(cl)
                xs, ys = coordinates(ln)
                # println(length(xs))
                vs     = [ W_full[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]
                # faint guide line
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)

                plot!(plt, xs .- 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)


                # faint guide line
                plot!(plt, xs, ys; lc=:black, lw=1, alpha=0.15, label=false)
                # weight‑coloured overlay
                plot!(plt, xs, ys; lc=1, lw=2, linealpha=vs, label=false)

                plot!(plt, xs, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                
                plot!(plt, xs, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                

                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)

                plot!(plt, xs .+ 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
            end
        end
    end

    plot!(plt, [-pi; -pi], [pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [-pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi; pi], [-pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi; 0], [pi; 0], lc = :black, lw = 1, label = false, ls = :dashdot)
    scatter!(plt, [π], [π], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [π], [0], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [0], [π], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [0], [0], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    # display(plt)
    return plt
end

####################################################################
############################## SOC + H #############################
####################################################################

" This function computes the DOS using the provided data obtained directly from the PROCAR file."

function plot_FS_ncollinear_H(kstart::Int=200, kend::Int=300, ΔE_F::Float64=0.0)
    # ------------------------------------------------------------
    # 1.  Read the slice we need
    # ------------------------------------------------------------
    kslice = kstart:kend
    h5f    = h5open(joinpath(FILE_DATA_SOC_H, "FS_H_201x201.h5"),"r")
    E_raw  = read(h5f,"eigenvalues")[:,:,kslice]
    E_F    = 5.53348380 - ΔE_F
    Energies = E_raw .- E_F

    layer_projs = sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4)
    Hlayer_proj = read(h5f,"proj_H_up")[:,:,kslice] .+ 
                  read(h5f,"proj_H_dn")[:,:,kslice]

    close(h5f)

    # sizes
    nkx, nky, nbands = size(Energies)

    # ------------------------------------------------------------
    # 2.  Mirror the first quadrant → full BZ helper
    # ------------------------------------------------------------
    mirror_q1(arr) = [reverse(arr)         reverse(arr; dims=1);
                      reverse(arr; dims=2) arr]

    E_full = cat((mirror_q1(@view Energies[:,:,b]) for b in 1:nbands)...; dims=3)
    W_full = cat((mirror_q1(@view layer_projs[:,:,b]) for b in 1:nbands)...; dims=3)

    H_full = cat((mirror_q1(@view Hlayer_proj[:,:,b]) for b in 1:nbands)...; dims=3) .* 15

    # k‑grid
    kx = range(-π, π; length=size(E_full,1))
    ky = range(-π, π; length=size(E_full,2))
    Δk = step(kx)

    # ------------------------------------------------------------
    # 3.  Small helper to map a coordinate → nearest grid index
    # ------------------------------------------------------------
    idx(x, g) = clamp(Int(round((x - first(g)) / Δk)) + 1, 1, length(g))

    # ------------------------------------------------------------
    # 4.  Build the plot
    # ------------------------------------------------------------
    plt = plot(; xlabel=L"\frac{k_x}{a}", ylabel = L"\frac{k_y}{a}",
                 xlims = (-6*π/5,6*π/5), ylims = (-6*π/5,6*π/5),
                 xticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                 yticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                framestyle=:box, aspect_ratio=:equal,
                 grid=false, size=(600,600))

    for b in 1:nbands
        for cl in levels(contours(kx,ky,E_full[:,:,b],[-0.00]))
            for ln in lines(cl)
                xs, ys = coordinates(ln)
                # println(length(xs))
                vs     = [ W_full[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]
                vs_H   = [ H_full[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]

                # faint guide line
                # plot!(plt, xs .- 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)

                # plot!(plt, xs .- 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- 2*π, ys; lc=2, lw=2, linealpha=vs_H, label=false)

                # plot!(plt, xs .- 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)

                # faint guide line
                # plot!(plt, xs, ys; lc=:black, lw=1, alpha=0.15, label=false)
                # weight‑coloured overlay
                plot!(plt, xs, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys; lc=2, lw=2, linealpha=vs_H, label=false)

                # plot!(plt, xs, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)
                
                # plot!(plt, xs, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)


                # plot!(plt, xs .+ 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)

                # plot!(plt, xs .+ 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=2, lw=2, linealpha=vs_H, label=false)
                
                # plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)
            end
        end
    end

    plot!(plt, [-pi; -pi], [pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [-pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi; pi], [-pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi; 0], [pi; 0], lc = :black, lw = 1, label = false, ls = :dashdot)
    scatter!(plt, [π], [π], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [π], [0], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [0], [π], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [0], [0], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    # display(plt)
    return plt
end

####################################################################
########################### SOC + H S217 ###########################
####################################################################

" This function computes the DOS using the provided data obtained directly from the PROCAR file."

function plot_FS_ncollinear_H_S(kstart::Int=400, kend::Int=600, ΔE_F::Float64=0.0)
    # ------------------------------------------------------------
    # 1.  Read the slice we need
    # ------------------------------------------------------------
    kslice = kstart:kend
    h5f    = h5open(joinpath(FILE_DATA_SOC_H_S, "FS_H_101x201.h5"),"r")
    E_raw  = read(h5f,"eigenvalues")[:,:,kslice]
    E_F    = 5.74715706 - ΔE_F
    Energies = E_raw .- E_F

    layer_projs = sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4)
    Hlayer_proj = read(h5f,"proj_H_up")[:,:,kslice] .+ 
                  read(h5f,"proj_H_dn")[:,:,kslice]

    close(h5f)

    # sizes
    nkx, nky, nbands = size(Energies)

    # ------------------------------------------------------------
    # 2.  Mirror the first quadrant → full BZ helper
    # ------------------------------------------------------------
    mirror_q1(arr) = [reverse(arr)         reverse(arr; dims=1);
                      reverse(arr; dims=2) arr]

    E_full = cat((mirror_q1(@view Energies[:,:,b]) for b in 1:nbands)...; dims=3)
    W_full = cat((mirror_q1(@view layer_projs[:,:,b]) for b in 1:nbands)...; dims=3)

    H_full = cat((mirror_q1(@view Hlayer_proj[:,:,b]) for b in 1:nbands)...; dims=3) .* 50

    # k‑grid
    kx = range(-π/2, π/2; length=size(E_full,1))
    ky = range(-π, π; length=size(E_full,2))
    Δk = step(kx)

    # ------------------------------------------------------------
    # 3.  Small helper to map a coordinate → nearest grid index
    # ------------------------------------------------------------
    idx(x, g) = clamp(Int(round((x - first(g)) / Δk)) + 1, 1, length(g))

    # ------------------------------------------------------------
    # 4.  Build the plot
    # ------------------------------------------------------------
    plt = plot(; xlabel=L"\frac{k_x}{a}", ylabel = L"\frac{k_y}{a}",
                 xlims = (-3*π/5,3*π/5), ylims = (-6*π/5,6*π/5),
                 xticks = ([-π/2:π/2:π/2;], ["-\\pi/2", "0", "\\pi/2"]),
                 yticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                framestyle=:box, aspect_ratio=:equal,
                 grid=false, size=(600,600))

    for b in 1:nbands
        for cl in levels(contours(kx,ky,E_full[:,:,b],[-0.00]))
            for ln in lines(cl)
                xs, ys = coordinates(ln)
                # println(length(xs))
                vs     = [ W_full[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]
                vs_H   = [ H_full[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]

                # faint guide line
                plot!(plt, xs .- π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- π, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)

                plot!(plt, xs .- π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- π, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- π, ys; lc=2, lw=2, linealpha=vs_H, label=false)

                plot!(plt, xs .- π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .- π, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)

                # faint guide line
                plot!(plt, xs, ys; lc=:black, lw=1, alpha=0.15, label=false)
                # weight‑coloured overlay
                plot!(plt, xs, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys; lc=2, lw=2, linealpha=vs_H, label=false)

                plot!(plt, xs, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)
                
                plot!(plt, xs, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)


                plot!(plt, xs .+ π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ π, ys .- 2*π; lc=2, lw=2, linealpha=vs_H, label=false)

                plot!(plt, xs .+ π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ π, ys; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ π, ys; lc=2, lw=2, linealpha=vs_H, label=false)
                
                plot!(plt, xs .+ π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
                plot!(plt, xs .+ π, ys .+ 2*π; lc=2, lw=2, linealpha=vs_H, label=false)
            end
        end
    end

    plot!(plt, [-pi/2; -pi/2], [pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi/2; pi/2], [-pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi/2; pi/2], [-pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi/2; pi/2], [pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi/2; 0], [pi; 0], lc = :black, lw = 1, label = false, ls = :dashdot)
    scatter!(plt, [π/2], [π], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [π/2], [0], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [0], [π], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [0], [0], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    # display(plt)
    return plt
end

function tst()
    println("Loading Fermi Surface Plots...")
end





####################################################################
############################## SOC #################################
####################################################################

" Test function to plot the Fermi Surface without weights. This is more similiar to an surface spectral function plot or a gapped map at E_F in the BZ"

function test_FS(kstart::Int=200, kend::Int=300, ΔE_F::Float64=0.0)
    # ------------------------------------------------------------
    # 1.  Read the slice we need
    # ------------------------------------------------------------
    kslice = kstart:kend
    h5f    = h5open(joinpath(FILE_DATA_SOC, "FS_201x201.h5"),"r")
    E_raw  = read(h5f,"eigenvalues")[:,:,kslice]
    E_F    = 5.62997313 - ΔE_F
    Energies = E_raw .- E_F

    layer_projs = sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4)
    # layer_projs = read(h5f, "proj1_dxy_up")[:,:,kslice] .+ read(h5f, "proj1_dxy_dn")[:,:,kslice]

    close(h5f)

    # sizes
    nkx, nky, nbands = size(Energies)
    min_ind = findmin(abs.(Energies), dims=3)
    E_tmp = Energies[min_ind[2]][:,:,1]
    

    # ------------------------------------------------------------
    # 2.  Mirror the first quadrant → full BZ helper
    # ------------------------------------------------------------
    mirror_q1(arr) = [reverse(arr)         reverse(arr; dims=1);
                      reverse(arr; dims=2) arr]

    E_full = mirror_q1(E_tmp)'

    color_map = :cork;
    # k‑grid
    kx = range(-π, π; length=size(E_full,1))
    ky = range(-π, π; length=size(E_full,2))
    Δk = step(kx)

    plt = heatmap(kx, ky, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  xlabel = L"\frac{k_x}{a}", ylabel = L"\frac{k_y}{a}",
                  xlims = (-6*π/5,6*π/5), ylims = (-6*π/5,6*π/5),
                  xticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                  yticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                  framestyle=:box, aspect_ratio=:equal,
                  grid=false, size=(600,600)
                  )
    plt = heatmap!(kx .- 2*π, ky, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    plt = heatmap!(kx .+ 2*π, ky, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )

    plt = heatmap!(kx, ky .- 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    plt = heatmap!(kx, ky .+ 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )

    plt = heatmap!(kx .- 2*π, ky .- 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    plt = heatmap!(kx .+ 2*π, ky .- 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    plt = heatmap!(kx .- 2*π, ky .+ 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    plt = heatmap!(kx .+ 2*π, ky .+ 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    
    # ------------------------------------------------------------
    kslice = kstart:kend
    h5f    = h5open(joinpath(FILE_DATA_SOC, "FS_201x201.h5"),"r")
    E_raw  = read(h5f,"eigenvalues")[:,:,kslice]
    E_F    = 5.62997313 - ΔE_F
    Energies = E_raw .- E_F

    layer_projs = sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4)
    # layer_projs = read(h5f, "proj1_dxy_up")[:,:,kslice] .+ read(h5f, "proj1_dxy_dn")[:,:,kslice]

    close(h5f)

    E_full = cat((mirror_q1(@view Energies[:,:,b]) for b in 1:nbands)...; dims=3)
    W_full = cat((mirror_q1(@view layer_projs[:,:,b]) for b in 1:nbands)...; dims=3)

    # ------------------------------------------------------------
    # 3.  Small helper to map a coordinate → nearest grid index
    # ------------------------------------------------------------
    idx(x, g) = clamp(Int(round((x - first(g)) / Δk)) + 1, 1, length(g))

    # ------------------------------------------------------------
    # 4.  Build the plot
    # ------------------------------------------------------------
    for b in 1:nbands
        for cl in levels(contours(kx,ky,E_full[:,:,b],[-0.00]))
            for ln in lines(cl)
                xs, ys = coordinates(ln)
                # println(length(xs))
                vs     = [ W_full[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]

                # faint guide line
                # plot!(plt, xs .- 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)

                # plot!(plt, xs .- 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)

                # plot!(plt, xs .- 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)

                # faint guide line
                # plot!(plt, xs, ys; lc=:black, lw=1, alpha=0.15, label=false)
                # weight‑coloured overlay
                plot!(plt, xs, ys; lc=1, lw=2, linealpha=vs, label=false)

                # plot!(plt, xs, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                
                # plot!(plt, xs, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)


                # plot!(plt, xs .+ 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)

                # plot!(plt, xs .+ 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                
                # plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
            end
        end
    end
    plot!(plt, [-pi; -pi], [pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [-pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi; pi], [-pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi; 0], [pi; 0], lc = :black, lw = 1, label = false, ls = :dashdot)
    scatter!(plt, [π], [π], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [π], [0], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [0], [π], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [0], [0], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    return plt
end

function test_FS_COHP(kstart::Int=200, kend::Int=300, ΔE_F::Float64=0.0)
    # ------------------------------------------------------------
    # 1.  Read the slice we need
    # ------------------------------------------------------------
    kslice = kstart:kend
    h5f    = h5open(joinpath(FILE_DATA_SOC, "FS_201x201.h5"),"r")
    E_raw  = read(h5f,"eigenvalues")[:,:,kslice]
    E_F    = 5.62997313 - ΔE_F
    Energies = E_raw .- E_F

    layer_projs = sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4)
    # layer_projs = read(h5f, "proj1_dxy_up")[:,:,kslice] .+ read(h5f, "proj1_dxy_dn")[:,:,kslice]

    close(h5f)

    # sizes
    nkx, nky, nbands = size(Energies)
    min_ind = findmin(abs.(Energies), dims=3)
    E_tmp = Energies[min_ind[2]][:,:,1]
    

    # ------------------------------------------------------------
    # 2.  Mirror the first quadrant → full BZ helper
    # ------------------------------------------------------------
    mirror_q1(arr) = [reverse(arr)         reverse(arr; dims=1);
                      reverse(arr; dims=2) arr]

    E_full = mirror_q1(E_tmp)'

    color_map = :cork;
    # k‑grid
    kx = range(-π, π; length=size(E_full,1))
    ky = range(-π, π; length=size(E_full,2))
    Δk = step(kx)

    plt = heatmap(kx, ky, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  xlabel = L"\frac{k_x}{a}", ylabel = L"\frac{k_y}{a}",
                  xlims = (-6*π/5,6*π/5), ylims = (-6*π/5,6*π/5),
                  xticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                  yticks = ([-π:π:π;], ["-\\pi", "0", "\\pi"]),
                  framestyle=:box, aspect_ratio=:equal,
                  grid=false, size=(600,600)
                  )
    plt = heatmap!(kx .- 2*π, ky, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    plt = heatmap!(kx .+ 2*π, ky, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )

    plt = heatmap!(kx, ky .- 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    plt = heatmap!(kx, ky .+ 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )

    plt = heatmap!(kx .- 2*π, ky .- 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    plt = heatmap!(kx .+ 2*π, ky .- 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    plt = heatmap!(kx .- 2*π, ky .+ 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    plt = heatmap!(kx .+ 2*π, ky .+ 2*π, E_full, lw = 0.0,
                  fill=true, levels = 1000,
                  color =  color_map, #:balance,:berlin,
                  )
    
    # ------------------------------------------------------------
    kslice = kstart:kend
    h5f    = h5open(joinpath(FILE_DATA_SOC, "FS_201x201.h5"),"r")
    E_raw  = read(h5f,"eigenvalues")[:,:,kslice]
    E_F    = 5.62997313 - ΔE_F
    Energies = E_raw .- E_F

    layer_projs = sum(read(h5f, "proj$(i)")[:,:,kslice] for i in 1:4)
    # layer_projs = read(h5f, "proj1_dxy_up")[:,:,kslice] .+ read(h5f, "proj1_dxy_dn")[:,:,kslice]

    close(h5f)

    E_full = cat((mirror_q1(@view Energies[:,:,b]) for b in 1:nbands)...; dims=3)
    W_full = cat((mirror_q1(@view layer_projs[:,:,b]) for b in 1:nbands)...; dims=3)

    # ------------------------------------------------------------
    # 3.  Small helper to map a coordinate → nearest grid index
    # ------------------------------------------------------------
    idx(x, g) = clamp(Int(round((x - first(g)) / Δk)) + 1, 1, length(g))

    # ------------------------------------------------------------
    # 4.  Build the plot
    # ------------------------------------------------------------
    for b in 1:nbands
        for cl in levels(contours(kx,ky,E_full[:,:,b],[-0.00]))
            for ln in lines(cl)
                xs, ys = coordinates(ln)
                # println(length(xs))
                vs     = [ W_full[idx(x,kx), idx(y,ky), b] for (x,y) in zip(xs,ys) ]

                # faint guide line
                # plot!(plt, xs .- 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)

                # plot!(plt, xs .- 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)

                # plot!(plt, xs .- 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .- 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)

                # faint guide line
                # plot!(plt, xs, ys; lc=:black, lw=1, alpha=0.15, label=false)
                # weight‑coloured overlay
                plot!(plt, xs, ys; lc=1, lw=2, linealpha=vs, label=false)

                # plot!(plt, xs, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)
                
                # plot!(plt, xs, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)


                # plot!(plt, xs .+ 2*π, ys .- 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .- 2*π; lc=1, lw=2, linealpha=vs, label=false)

                # plot!(plt, xs .+ 2*π, ys; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys; lc=1, lw=2, linealpha=vs, label=false)
                
                # plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=:black, lw=1, alpha=0.15, label=false)
                plot!(plt, xs .+ 2*π, ys .+ 2*π; lc=1, lw=2, linealpha=vs, label=false)
            end
        end
    end
    plot!(plt, [-pi; -pi], [pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [-pi; -pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi; pi], [-pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [-pi; pi], [pi; pi], lc = :black, lw = 2, label = false)
    plot!(plt, [pi; 0], [pi; 0], lc = :black, lw = 1, label = false, ls = :dashdot)
    scatter!(plt, [π], [π], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [π], [0], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [0], [π], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    scatter!(plt, [0], [0], label = false, ms = 5, markerstrokewidth = 0.2, color = 3)
    return plt
end