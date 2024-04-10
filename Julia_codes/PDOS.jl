using Plots
using DelimitedFiles
using CSV

const FOLDER = raw"/home/r_floren/BandStructure/QE/NbP/Main_folder/supercell/pdos/supercell/"

function openfile(data::String)
    datafile = readdlm(data)
    energy = datafile[2:end,1]
    pdos = datafile[2:end,3]
    return energy, pdos
end

const FILE_FE = raw"/home/r_floren/BandStructure/QE/NbP/Main_folder/supercell/pdos/supercell/NbP.supercell2x2.nscf.out"

function openfile1(data::String)
    open(data, "r") do file
        return readlines(file)
    end
end

function extract_fermi_energy(FILE_PATH::String)::Float64
    lines = openfile1(FILE_PATH)
    match = findfirst(i -> occursin("Fermi", i), lines)
    if match isa Nothing
        error("Fermi energy not found in the input data.")
    else
        return parse(Float64, split(lines[match])[end-1])
    end
end

# const files = filter(x->occursin("NbP\\.supercell2x2\\.pdos\\.dat\\.pdos_atm#\\d(Nb)_wfc#.(d)",x), readdir("/home/r_floren/BandStructure/QE/NbP/Main_folder/supercell/pdos/supercell/"))
const files = cat(filter(x->occursin(r"NbP\.supercell2x2\.pdos\.dat\.pdos_atm#.\(Nb\)_wfc#.\(d\)",x), readdir(FOLDER)),
filter(x->occursin(r"NbP\.supercell2x2\.pdos\.dat\.pdos_atm#.\(P\)_wfc#.\(p\)",x), readdir(FOLDER)),
filter(x->occursin(r"NbP\.supercell2x2\.pdos\.dat\.pdos_atm#\d{2}\(P\)_wfc#.\(p\)",x), readdir(FOLDER)),dims=1)


function plot_pdos(FOLDER_FILES::String,FILES::Vector{String},M::Vector{Int64})
    p = plot()
    clrs = [:blue, :red, :green, :purple, :orange];
    linestyle = [:solid, :dash, :dashdot, :dot]
    for i in FILES[[1,6,9,13]]
        E,  PDOS = openfile(FOLDER_FILES*i)
        region_condition = zeros(length(E))
        region_condition[E .< extract_fermi_energy(FILE_FE)] = PDOS[E .< extract_fermi_energy(FILE_FE)]

        plot!(p,E, PDOS,
        lw=1.5, label=i[36:end-9]*"-"*i[end-3:end],
        lc=clrs[findall(x->x==i,FILES[[1,6,9,13]])],
        ls=linestyle[findall(x->x==i,FILES[[1,6,9,13]])])
        # plot!(p,E, PDOS,
        # lw=1.5, label=i[36:end-9]*"-"*i[end-3:end], 
        # ls = linestyle[parse(Int, string(i[end-3]))],
        # lc=clrs[parse(Int, string(i[end-3]))])

        plot!(p,E, region_condition,
        fill=(0, region_condition,clrs[findall(x->x==i,FILES[[1,6,9,13]])], 0.15),
        label=false,lw=0)

        # plot!(p,E, region_condition,
        # fill=(0, region_condition,clrs[parse(Int, string(i[end-3]))], 0.15),
        # lw=0, label=false,
        # ls = linestyle[parse(Int, string(i[36]))])
    end
    Enn,  PPDOS = openfile(FOLDER_FILES*filter(x->occursin(r".pdos_tot",x), readdir(FOLDER_FILES))[1])
    plot!(p,extract_fermi_energy(FILE_FE)*ones(2), [0, findmax(PPDOS)[1]],lw=0.75,lc=:black,ls=:dashdot,label=false)
    plot!(p,extract_fermi_energy(FILE_FE)*ones(2), [0, 2],lw=0.75,lc=:black,ls=:dashdot,label=false)
    plot!(p,grid=false,
    yticks=false, xlabel="Energy (eV)",
    ylabel="DOS (states/eV)", xlims=(M[1],M[2]),
    ylims=(0, 2), legend=:topleft)
    #ylims=(0, findmax(PPDOS)[1]-10), legend=:topleft)
    annotate!(p,extract_fermi_energy(FILE_FE)-0.1,1.5,
    Plots.text("Fermi energy", 12, :dark, rotation = 90 ) )
    plot!(twinx(),Enn, PPDOS, lw=2,label="Total DOS", lc=:black,
    xlims=(M[1],M[2]), ylims=(0, 50),ylabel="Total DOS (states/eV)")

    return p
end
