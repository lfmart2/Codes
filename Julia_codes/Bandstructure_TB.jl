using Plots
using DelimitedFiles
using HDF5
using CSV
using Plots.PlotMeasures
using DataFrames
using LinearAlgebra
using LaTeXStrings

## Files for the VASP  data
const PATH_DF_nSpin_nH         = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Bandstructure/df/"
const FILE_BANDS_INFO_nSpin_nH = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Bandstructure/df/KPOINTS.tot"

const PATH_DF_nSpin_H          = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_H/Bandstructure/df/"
const FILE_BANDS_INFO_nSpin_H  = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_H/Bandstructure/df/KPOINTS.tot"

const PATH_DF_Spin_nH          = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/df/"
const FILE_BANDS_INFO_Spin_nH  = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/Bandstructure/df/KPOINTS.tot"

const PATH_DF_Spin_H           = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/df/"
const FILE_BANDS_INFO_Spin_H   = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/Bandstructure/df/KPOINTS.tot"

## Files For the Tight-Binding data. Bandstructure calculatated from Wannier90
const FILE_BANDS_TB_nSpin_nH   = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/wannier90/NbP_band.dat"
const FILE_TB_INFO_nSpin_nH    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/wannier90/NbP_band.labelinfo.dat"

## Files For the Tight-Binding data. Bandstructure calculatated from Wannier90
const FILE_BANDS_TB_nSpin_nH1   = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/wannier90_2/NbP_band.dat"
const FILE_TB_INFO_nSpin_nH1    = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/Tight-Binding/wannier90_2/NbP_band.labelinfo.dat"
# Pre Proccesing functions for the SCF (VASP) and Tight-Bidning Bandstructure:

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

function label_sym_points(FILE_PATH::String)
    data = readdlm(FILE_PATH)[5:end,end]
    x_label = String[]
    if data[1] == "Gamma"
        x_label = push!(x_label,"\$ \\"* data[1]* "\$")
    else
        x_label = push!(x_label,"\$"* data[1]* "\$")
    end

    for i=2:2:length(data)-1
        if data[i] == data[i+1]
            if data[i] == "Gamma"
                x_label = push!(x_label,"\$ \\"* data[i]* "\$")
            else
                x_label = push!(x_label,"\$"* data[i]* "\$")
            end
        else
            if data[i] == "Gamma"
                x_label = push!(x_label,"\$ \\"* data[i]* "|"* data[i+1]* "\$")
            elseif data[i+1] == "Gamma"
                x_label = push!(x_label,"\$"* data[i]* "|\\"* data[i+1]* "\$")
            else
                x_label = push!(x_label,"\$"* data[i]* "|"* data[i+1]* "\$")
            end
        end
    end
    if data[end] == "Gamma"
        x_label = push!(x_label,"\$ \\"* data[end]* "\$")
    else
        x_label = push!(x_label,"\$"* data[end]* "\$")
    end

end


function label_sym_points_tb(FILE_PATH::String)
    data = readdlm(FILE_PATH)[1:end,1]
    x_label = String[]
    if data[1] == "Gamma" || data[1] == "G"
        x_label = push!(x_label,"\$ \\Gamma \$")
    else
        x_label = push!(x_label,"\$"* data[1]* "\$")
    end

    for i=2:1:length(data)-1
        if data[i] == "Gamma" || data[i] == "G"
            x_label = push!(x_label,"\$ \\Gamma \$")
        else
            x_label = push!(x_label,"\$"* data[i]* "\$")
        end
    end
    if data[end] == "Gamma" || data[end] == "G"
        x_label = push!(x_label,"\$ \\Gamma \$")
    else
        x_label = push!(x_label,"\$"* data[end]* "\$")
    end

end

function openfile(data::String)
    open(data, "r") do file
        return readlines(file)
    end
end


# Function to plot the bandstructure from the tight-binding calculation
function plot_band_tb(FILE_PATH_BAND_TB::String)
    # Range of Energy
    label = ["Y" 0; "\$ \\Gamma \$"  0.4679950372;"X"  0.9359900745;"M"  1.4039851117;"\$ \\Gamma \$"  2.0658300405]
    label[:,2] = label[:,2]./maximum(label[:,2])
    # Parameters for the plot
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:red, :blue, :black]

    fermi_en = 10.3786
    
    df_tb = open(FILE_PATH_BAND_TB) do io
        dropmissing(CSV.read(io, DataFrame; header=["k", "En"], ignoreemptyrows=true, delim=' ', ignorerepeated=true))
    end
    df_tb.k = df_tb.k ./ maximum(df_tb.k)

    p = plot(leftmargin = 1mm)    


    n_kpoints_tb =  findall(x->df_tb.k[x] == 0.0,1:length(df_tb.k))[1:2]
    n_kpoints_tb = n_kpoints_tb[2]-n_kpoints_tb[1]
    n_bnds_tb = length(df_tb.k)÷n_kpoints_tb

    for i in 1:n_bnds_tb
        plot!(p, df_tb.k[(i-1)*n_kpoints_tb+1:i*n_kpoints_tb], df_tb.En[(i-1)*n_kpoints_tb+1:i*n_kpoints_tb].-fermi_en,
            lw=LINE_WIDTH[1],
            lc=LINE_COLOR[2],
            label=false
            )
    end

    plot!(p, grid=false, xticks=false,
        ylabel=L"Energy - $E_F$  (eV)",
        ylims=(-1.0, 1.0),
        xlims=(0.0,findmax(df_tb.k)[1]),
        )

    hline!(p, [0],
        lw=0.75,
        lc=LINE_COLOR[3],
        label=false
        )
    vline!(p, label[:,2],
        lw=LINE_WIDTH[4],
        lc=:black,
        ls=:solid,label=false
        )
        plot!(p,xticks = (label[:,2],label[:,1]),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 16,
        legendfontsize = 12)

    plot!(p, df_tb.k[1:n_kpoints_tb], df_tb.En[1:n_kpoints_tb].-fermi_en,
            lw=LINE_WIDTH[1],
            lc=LINE_COLOR[2],
            ls=:dashdot,
            label="Tight-binding"
        )
    return p
end


# Function to plot the projected bandstructure from the SCF (VASP) calculation. The output are the Parameters
# for a plot with the bandstructure from the SCF and the tight-binding calculation.
function proj_nBandstructure_vasp(FILE_PATH_BS::String,FILE_PATH_SYM::String)
    p1 = plot();
     # Parameters for the plot
     LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
     LINE_COLOR = [:red, :blue, :black]

 
    fermi_energy = CSV.File(open(FILE_PATH_BS*"df_energies_total.csv"); limit=5).energy[5]
    tot_atoms = CSV.File(open(FILE_PATH_BS*"df_energies_total.csv"); limit=4).energy[4]

    sym_points = readdlm(FILE_PATH_SYM)[5:end,1:3]

    #Calculation of the energies
    n_kpoints, n_bnds, energies_2, data_2 = load_filtered_data(FILE_PATH_BS*"df_bands_total.csv", FILE_PATH_BS*"df_energies_total.csv", ["s","p","d"], [3,6,9,12,51,54,57,60])
    tot_atoms = 8; # Number of atoms in the unit cell. Delete for Total Bandstructure
    tmp = sum(eachcol(data_2))
    tmp_reshape = reshape(tmp, tot_atoms,:)
    tmp = vec(sum(tmp_reshape, dims=1))
    bands_2 = reshape(tmp, n_bnds, n_kpoints)'

    # Calculate sections and x_axis
    secs = cumsum([0;norm.(eachrow(diff(sym_points,dims=1)),2)[1:2:end]])
    secs = secs ./ maximum(secs)
    x_axis = vcat(range.(secs[1:end-2],secs[2:end-1],length=50)...,range(secs[end-1],secs[end], length = 70)...)
    x_axis = x_axis ./ maximum(x_axis)
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
    # plot!(p1,grid=false, ylabel="Energy (eV)", ylim=(-2.5,3.5))
    plot!(p1,grid=false, ylabel=L"Energy - $E_F$  (eV)",
        ylims=(-1,1),
        xlims=(0,x_axis[end]),
    )
    plot!(p1, xticks = (secs,label_sym_points(FILE_PATH_SYM)),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 16,
        legendfontsize = 12)
    display(p1)
    return x_axis, reshape(energies_2, n_bnds, n_kpoints)'.-fermi_energy, bands_2, secs
end

# Function to plot the projected bandstructure from the tight-binding calculation. The output are the Parameters
# for a plot with the bandstructure from the SCF and the tight-binding calculation.
function plot_proj_band_tb(FILE_PATH_BAND_TB::String,FILE_PATH_INFO_TB::String)
    # Range of Energy
    label = label_sym_points_tb(FILE_PATH_INFO_TB)
    label_x = readdlm(FILE_PATH_INFO_TB)[1:end,3]
    label_x = label_x./maximum(label_x)
    # Parameters for the plot
    LINE_WIDTH = [2.0, 2.0, 0.75, 1.0]
    LINE_COLOR = [:red, :blue, :black]

    fermi_en = 10.3786
    
    df_tb = open(FILE_PATH_BAND_TB) do io
        dropmissing(CSV.read(io, DataFrame; header=["k", "En", "proj"], ignoreemptyrows=true, delim=' ', ignorerepeated=true))
    end
    df_tb.k = df_tb.k ./ maximum(df_tb.k)

    p = plot(leftmargin = 1mm)    


    n_kpoints_tb =  findall(x->df_tb.k[x] == 0.0,1:length(df_tb.k))[1:2]
    n_kpoints_tb = n_kpoints_tb[2]-n_kpoints_tb[1]
    n_bnds_tb = length(df_tb.k)÷n_kpoints_tb

    for i in 1:n_bnds_tb
        plot!(p, df_tb.k[(i-1)*n_kpoints_tb+1:i*n_kpoints_tb], df_tb.En[(i-1)*n_kpoints_tb+1:i*n_kpoints_tb].-fermi_en,
            lw=LINE_WIDTH[1],
            alpha = df_tb.proj[(i-1)*n_kpoints_tb+1:i*n_kpoints_tb],
            lc=LINE_COLOR[2],
            label=false
            )
    end

    plot!(p, grid=false, xticks=false,
        ylabel=L"Energy - $E_F$  (eV)",
        ylims=(-1.0, 1.0),
        xlims=(0.0,findmax(df_tb.k)[1]),
        )

    hline!(p, [0],
        lw=LINE_WIDTH[3],
        lc=LINE_COLOR[3],
        label=false,
        ls = :dashdot
        )
    vline!(p, label_x,
        lw=LINE_WIDTH[4],
        lc=:black,
        ls=:solid,label=false
        )
        plot!(p,xticks = (label_x,label),
        ytickfontsize  = 10,
        xtickfontsize  = 12,
        yguidefontsize = 16,
        legendfontsize = 12)

    plot!(p, df_tb.k[1:n_kpoints_tb], df_tb.En[1:n_kpoints_tb].-fermi_en,
            lw=LINE_WIDTH[1],
            lc=LINE_COLOR[2],
            ls=:dashdot,
            label="Tight-binding"
        )
    # display(p)
    return n_kpoints_tb, n_bnds_tb, df_tb.k, df_tb.En .- fermi_en, df_tb.proj, label
end

# Function to plot the bandstructure from the SCF (VASP) and the tight-binding calculation, merged in the same plot
function plot_tb_scf()
    # Parameters for the plot
    LINE_WIDTH = [2.25, 2.0, 0.75, 1.0]
    LINE_COLOR = [:red, :blue, :black]

    x_axis_scf, y_axis_scf, alpha_scf, secs_scf = proj_nBandstructure_vasp(PATH_DF_nSpin_nH,FILE_BANDS_INFO_nSpin_nH)
    n_kpoints_tb, n_bnds_tb, x_axis_tb, y_axis_tb, alpha_tb, label = plot_proj_band_tb(FILE_BANDS_TB_nSpin_nH1,FILE_TB_INFO_nSpin_nH1)
    p1 = plot()

    plot!(p1,x_axis_scf,y_axis_scf,
        markersize=LINE_WIDTH[1],
        seriescolor=LINE_COLOR[1],
        seriestype = :scatter,
        seriesalpha = alpha_scf,
        label=false,xticks=false
        )
    plot!(p1,(x_axis_scf[1],y_axis_scf[1]),
        markersize=LINE_WIDTH[1],
        seriescolor=LINE_COLOR[1],
        seriestype = :scatter,
        label="VASP"
        )
    
    # Plot fermi energy and sections
    hline!(p1,[0],
        lw=LINE_WIDTH[3],
        lc=LINE_COLOR[3],
        ls=:dashdot,label=false
        )
    vline!(p1,secs_scf,lw=LINE_WIDTH[4],
        lc=:black,
        ls=:solid,label=false
        )
    plot!(p1,grid=false, ylabel=L"Energy - $E_F$  (eV)",
        ylims=(-1,1),
        xlims=(0,x_axis_scf[end]),
    )
    plot!(p1, xticks = (secs_scf,label),
    ytickfontsize  = 10,
    xtickfontsize  = 12,
    yguidefontsize = 16,
    legendfontsize = 12)

    for i in 1:n_bnds_tb
        plot!(p1, x_axis_tb[(i-1)*n_kpoints_tb+1:i*n_kpoints_tb], y_axis_tb[(i-1)*n_kpoints_tb+1:i*n_kpoints_tb].-0.05,
            lw=LINE_WIDTH[1],
            alpha = alpha_tb[(i-1)*n_kpoints_tb+1:i*n_kpoints_tb],
            lc=LINE_COLOR[2],
            label=false
            )
    end
    # plot!(p1,range(0, 1, length = length(kx)),Energies' .- 0.05,
    #     lw=LINE_WIDTH[1],
    #     seriesalpha = Proj',
    #     lc=LINE_COLOR[2],
    #     label=false
    #     )
    plot!(p1, x_axis_tb[1:n_kpoints_tb], y_axis_tb[1:n_kpoints_tb],
    lw=LINE_WIDTH[1],
    lc=LINE_COLOR[2],
    ls=:dashdot,
    label="Tight-binding"
    )
    return p1
end