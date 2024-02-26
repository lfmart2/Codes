using Plots
using DelimitedFiles
using CSV



function openfile(data::String)
    datafile = readdlm(data)
    energy = datafile[2:end,1]
    pdos = datafile[2:end,3]
    return energy, pdos
end

const FILE_FE = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/Local/Bulk/scf.out"

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

const files = filter(x->occursin("NbP_pdos.dat.pdos_atm",x), readdir("/home/r_floren/BandStructure/QE/NbP/github/Outputs/fire_cluster/NbP_spin/"))


function plot_pdos(Folder_PATH::Vector{String},M::Vector{Int64})
    p = plot()
    clrs = [:blue, :red, :green, :blue];
    linestyle = [:solid, :dash, :dashdot, :dot]
    for i in Folder_PATH
        E,  PDOS = openfile("/home/r_floren/BandStructure/QE/NbP/github/Outputs/fire_cluster/NbP_spin/"*i)
        region_condition = zeros(length(E))
        region_condition[E .> extract_fermi_energy(FILE_FE)] = PDOS[E .> extract_fermi_energy(FILE_FE)]
        
        plot!(p,E, PDOS,
        lw=1.5, label=i[25:end-10]*" "*i[end-3:end], 
        ls = linestyle[parse(Int, string(i[23]))],
        lc=clrs[parse(Int, string(i[end-3]))])

        plot!(p,E, PDOS,
        fill=(0, region_condition,clrs[parse(Int, string(i[end-3]))], 0.25),
        lw=1.5, label=false, 
        ls = linestyle[parse(Int, string(i[23]))],
        lc=clrs[parse(Int, string(i[end-3]))])
    end
    Enn,  PPDOS = openfile("/home/r_floren/BandStructure/QE/NbP/github/Outputs/fire_cluster/NbP_spin/NbP_pdos.dat.pdos_tot")
    plot!(p,extract_fermi_energy(FILE_FE)*ones(2), [0, findmax(PPDOS)[1]],lw=0.75,lc=:black,ls=:dashdot,label=false)
    plot!(p,grid=false, yticks=false, xlabel="Energy (eV)", ylabel="DOS (states/eV)", xlims=(M[1],M[2]), ylims=(0, findmax(PPDOS)[1]-10), legend=:topleft)
    plot!(p,Enn, PPDOS, lw=2,label="Total DOS", lc=:black)
    annotate!(p,extract_fermi_energy(FILE_FE)-0.2,4, Plots.text("Fermi energy", 12, :dark, rotation = 90 ))
    return p
end
