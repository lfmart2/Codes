using Plots
using DelimitedFiles
using CSV
using Plots.PlotMeasures
using DataFrames
using LinearAlgebra
using LaTeXStrings


## Files to lead the QE output data
const FILE_DOS1 = raw"/home/r_floren/BandStructure/VASP/NbP/Supercell_final/Coord_1/nSpin_nH/DOS/DOSCAR"
const FILE_DOS2 = raw"/home/r_floren/BandStructure/VASP/NbP/Supercell_final/Coord_1/nSpin_H/DOS/DOSCAR"
const FILE_DOS3 = raw"/home/r_floren/BandStructure/VASP/NbP/Supercell_final/Coord_1/Spin_nH/DOS/DOSCAR"
const FILE_DOS4 = raw"/home/r_floren/BandStructure/VASP/NbP/Supercell_final/Coord_1/Spin_H/DOS/DOSCAR"


function openfile(data::String)
    open(data, "r") do file
        return readlines(file)
    end
end

function extract_fermi_energy(FILE_PATH::String)::Float64
    return parse(Float64, split(openfile(FILE_PATH)[6])[4])
end

function natomic_orbital(a_orb::Vector{String})
    a_orbital = String[]
    # a_orbital = push!(a_orbital,"En")
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

function spin_atomic_orbital(spin::String, a_orb::Vector{String})
    a_orbital = String[]
    # a_orbital = push!(a_orbital,"En")
    if spin == "up"
        for i in a_orb
            if i == "s"
                a_orbital = push!(a_orbital,"u-s")
            elseif i == "p"
                a_orbital = push!(a_orbital,"u-px","u-py","u-pz")
            elseif i == "d"
                a_orbital = push!(a_orbital,"u-dxy","u-dyz","u-dz2","u-dxz","u-dx2")
            end
        end
    elseif spin == "down"
        for i in a_orb
            if i == "s"
                a_orbital = push!(a_orbital,"d-s")
            elseif i == "p"
                a_orbital = push!(a_orbital,"d-px","d-py","d-pz")
            elseif i == "d"
                a_orbital = push!(a_orbital,"d-dxy","d-dyz","d-dz2","d-dxz","d-dx2")
            end
        end
    end
    return a_orbital
end


function plot_pdos(FILE_PATH_DOSCAR::String,a_orbital::Vector{Vector{String}},atoms::Vector{Vector{Int64}})
    # Range of Energy
    #E_min = readdlm(FILE_PATH_DOSCAR)[6,2];
    #E_max = readdlm(FILE_PATH_DOSCAR)[6,1];
    # Parameters for the plot
    LINE_WIDTH = 1.5
    LINE_COLOR = [:red, :blue, :green, :black]
    ANNOTATE_POSITION = -0.05
    ANNOTATE_SIZE = 12
    fermi_energy = extract_fermi_energy(FILE_PATH_DOSCAR)
    println(fermi_energy)
    e_mesh = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[6])[3])
    n_atoms = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[1])[1])
    df = open(FILE_PATH_DOSCAR) do io
        dropmissing(CSV.read(io, DataFrame; header=["En", "s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2"], skipto =e_mesh+7, ignoreemptyrows=true, delim=' ', ignorerepeated=true),:dx2)
    end
    energy_data = df[4002:2*4001,1] .- fermi_energy

    PDOS = zeros(Float64, size(atoms)[1], e_mesh)


    for j in 1:size(atoms)[1]
        for i in atoms[j]
            PDOS[j,:] += sum(eachcol(select(df[(i-1)*e_mesh+1:i*e_mesh,:], natomic_orbital(a_orbital[j]) ) ))
        end
    end

    p = plot(leftmargin = 5mm)
    for j in 1:size(atoms)[1]
        plot!(p, PDOS[j,:], energy_data,
            lw=LINE_WIDTH,
            lc=j,
            label=a_orbital[j]
            )
    end

    plot!(p, grid=false, xticks=false,
        ylabel="Energy (eV)",
        xlabel="DOS (states/eV)",
        ylims=(-5.0, 4.0),
        xlims=(0.0,9.4652) #findmax(PDOS[1,50:end])[1]+1),
        # yticks = ([0],[L"$E_F$"]),
        )

    hline!(p, [0],
        lw=0.75,
        lc=LINE_COLOR[4],
        ls=:dashdot,
        label=false
        )
    # plot!(p, xlims=(fermi_energy-2, fermi_energy+2),ylims=(0, 7))
    annotate!(p, ANNOTATE_POSITION, 15, Plots.text("Fermi energy", ANNOTATE_SIZE, :dark, rotation = 90 ))
    return p
end


function plot_pdos2(FILE_PATH_DOSCAR::String,a_orbital::Vector{Vector{String}},atoms::Vector{Vector{Int64}})
    # Range of Energy
    #E_min = readdlm(FILE_PATH_DOSCAR)[6,2];
    #E_max = readdlm(FILE_PATH_DOSCAR)[6,1];
    # Parameters for the plot
    LINE_WIDTH = 1.5
    LINE_COLOR = [:black, :gray, :green, :black]
    ANNOTATE_POSITION = -0.05
    ANNOTATE_SIZE = 12
    LABEL = ["d-orbitals","s-orbitals","Hydrogen"]
    fermi_energy = extract_fermi_energy(FILE_PATH_DOSCAR)
    println(fermi_energy)
    e_mesh = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[6])[3])
    n_atoms = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[1])[1])
    df = open(FILE_PATH_DOSCAR) do io
        dropmissing(CSV.read(io, DataFrame; header=["En", "s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2"], skipto =e_mesh+7, ignoreemptyrows=true, delim=' ', ignorerepeated=true),:dx2)
    end
    energy_data = df[4002:2*4001,1] .- fermi_energy

    PDOS = zeros(Float64, size(atoms)[1], e_mesh)


    for j in 1:size(atoms)[1]
        for i in atoms[j]
            PDOS[j,:] += sum(eachcol(select(df[(i-1)*e_mesh+1:i*e_mesh,:], natomic_orbital(a_orbital[j]) ) ))
        end
    end

    p = plot(leftmargin = 5mm)
    for j in 1:size(atoms)[1]
        plot!(p, energy_data, PDOS[j,:],
            lw=LINE_WIDTH,
            lc=LINE_COLOR[j],
            label=LABEL[j]
            )
    end

    plot!(p, grid=false, yticks=false,
        xlabel="Energy (eV)",
        ylabel="DOS (states/eV)",
        xlims=(-4.0, 4.0),
        ylims=(0.0,9.4652),
        title = "No spin + H" #findmax(PDOS[1,50:end])[1]+1),
        # yticks = ([0],[L"$E_F$"]),
        )

    vline!(p, [0],
        lw=0.75,
        lc=LINE_COLOR[4],
        ls=:dashdot,
        label=false
        )
    # plot!(p, xlims=(fermi_energy-2, fermi_energy+2),ylims=(0, 7))
    annotate!(p, ANNOTATE_POSITION, 15, Plots.text("Fermi energy", ANNOTATE_SIZE, :dark, rotation = 90 ))
    return p
end




function plot_spin_pdos(FILE_PATH_DOSCAR::String,a_orbital::Vector{Vector{String}},atoms::Vector{Vector{Int64}},spn::Vector{String})
    # Range of Energy
    #E_min = readdlm(FILE_PATH_DOSCAR)[6,2];
    #E_max = readdlm(FILE_PATH_DOSCAR)[6,1];
    # Parameters for the plot
    LINE_WIDTH = 1.75
    LINE_COLOR = [:red, :blue, :black]
    ANNOTATE_POSITION = -0.05
    ANNOTATE_SIZE = 12
    fermi_energy = extract_fermi_energy(FILE_PATH_DOSCAR)
    e_mesh = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[6])[3])
    n_atoms = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[1])[1])
    df = open(FILE_PATH_DOSCAR) do io
        dropmissing(CSV.read(io, DataFrame; header=["En", "u-s", "d-s", "u-py", "u-pz", "u-px", "d-py", "d-pz", "d-px", "u-dxy", "u-dyz", "u-dz2", "u-dxz", "u-dx2", "d-dxy", "d-dyz", "d-dz2", "d-dxz", "d-dx2"], skipto =e_mesh+7, ignoreemptyrows=true, delim=' ', ignorerepeated=true),"d-dx2")
    end
    energy_data = df[4002:2*4001,1] .- fermi_energy
    PDOS1 = zeros(Float64, e_mesh)
    PDOS2 = zeros(Float64, e_mesh)
    for i in atoms[1]
        PDOS1 += sum(eachcol(select(df[(i-1)*e_mesh+1:i*e_mesh,:], spin_atomic_orbital(spn[1], a_orbital[1])) ))
        PDOS1 += sum(eachcol(select(df[(i-1)*e_mesh+1:i*e_mesh,:], spin_atomic_orbital(spn[2], a_orbital[1])) ))
    end
    for i in atoms[2]
        PDOS2 += sum(eachcol(select(df[(i-1)*e_mesh+1:i*e_mesh,:], spin_atomic_orbital(spn[1], a_orbital[2])) ))
        PDOS2 += sum(eachcol(select(df[(i-1)*e_mesh+1:i*e_mesh,:], spin_atomic_orbital(spn[2], a_orbital[2])) ))
    end

    p = plot(leftmargin = 5mm)
    plot!(p, PDOS1, energy_data,
        lw=LINE_WIDTH,
        lc=3,
        label="Total d orbitals"
        )
    plot!(p, PDOS2, energy_data,
        lw=LINE_WIDTH,
        lc=4,
        label="Total s-p orbitals"
        )
    plot!(p, grid=false, xticks=false,
        ylabel = "", # ylabel="Energy (eV)",
        xlabel="DOS",
        ylims=(-1.0, 1.0),
        xlims=(0.0,20),
        # xlims=(0.0,findmax(PDOS1[50:end])[1]),
        # yticks = ([0],[L"$E_F$"]),
        )

    hline!(p, [0],
        lw=0.75,
        lc=LINE_COLOR[3],
        ls=:dashdot,
        label=false
        )
    plot!(p,
        yticks = false,
        xtickfontsize  = 10,
        yguidefontsize = 18,
        legend = false
        )
    # plot!(p, xlims=(fermi_energy-2, fermi_energy+2),ylims=(0, 7))
    annotate!(p, ANNOTATE_POSITION, 15, Plots.text("Fermi energy", ANNOTATE_SIZE, :dark, rotation = 90 ))
    return p
end

function plot_spin_pdos2(FILE_PATH_DOSCAR::String,a_orbital::Vector{Vector{String}},atoms::Vector{Vector{Int64}})
    # Range of Energy
    #E_min = readdlm(FILE_PATH_DOSCAR)[6,2];
    #E_max = readdlm(FILE_PATH_DOSCAR)[6,1];
    # Parameters for the plot
    LINE_WIDTH = 1.75
    LINE_COLOR = [:black, :gray, :green, :black]
    LABEL = ["d-orbitals up","d-orbitals down","sp-orbitals up","sp-orbitals down","Hydrogen up","Hydrogen down"]
    ANNOTATE_POSITION = -0.05
    ANNOTATE_SIZE = 12
    fermi_energy = extract_fermi_energy(FILE_PATH_DOSCAR)
    e_mesh = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[6])[3])
    n_atoms = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[1])[1])
    df = open(FILE_PATH_DOSCAR) do io
        dropmissing(CSV.read(io, DataFrame; header=["En", "u-s", "d-s", "u-py", "u-pz", "u-px", "d-py", "d-pz", "d-px", "u-dxy", "u-dyz", "u-dz2", "u-dxz", "u-dx2", "d-dxy", "d-dyz", "d-dz2", "d-dxz", "d-dx2"], skipto =e_mesh+7, ignoreemptyrows=true, delim=' ', ignorerepeated=true),"d-dx2")
    end
    energy_data = df[4002:2*4001,1] .- fermi_energy

    PDOS = zeros(Float64, 2*size(atoms)[1], e_mesh)


    for j in 1:size(atoms)[1]
        for i in atoms[j]
            PDOS[2*j-1,:] += sum(eachcol(select(df[(i-1)*e_mesh+1:i*e_mesh,:], spin_atomic_orbital("up", a_orbital[j])) ))
            PDOS[2*j,:] += sum(eachcol(select(df[(i-1)*e_mesh+1:i*e_mesh,:], spin_atomic_orbital("down", a_orbital[j])) ))
        end
    end

    p = plot(leftmargin = 5mm)
    for j in 1:size(atoms)[1]
        plot!(p, energy_data, PDOS[2*j-1,:],
            lw=LINE_WIDTH,
            ls = :solid,
            lc=LINE_COLOR[j],
            label=LABEL[2*j-1]
            )
        plot!(p, energy_data, PDOS[2*j-1,:],
            lw=LINE_WIDTH,
            ls = :dash,
            lc=LINE_COLOR[j],
            label=LABEL[2*j]
            )
    end
    # plot!(p, energy_data, PDOS11,
    #     lw=LINE_WIDTH,
    #     ls = :solid,
    #     lc=1,
    #     label="Spin-up d orbitals"
    #     )
    # plot!(p, energy_data, PDOS12, 
    #     lw=LINE_WIDTH,
    #     ls = :dash,
    #     lc=1,
    #     label="Spin-down d orbitals"
    #     )
    # plot!(p, energy_data, PDOS21,
    #     lw=LINE_WIDTH,
    #     ls = :solid,
    #     lc=2,
    #     label="Spin-up s-p orbitals"
    #     )
    # plot!(p, energy_data, PDOS22, 
    #     lw=LINE_WIDTH,
    #     ls = :dash,
    #     lc=2,
    #     label="Spin-down s-p orbitals"
    #     )
    plot!(p, grid=false, yticks=false,
        xlabel="Energy (eV)",
        ylabel="DOS (states/eV)",
        xlims=(-4.0, 4.0),
        ylims=(0.0,5),
        title = "Spin + H"
        # xlims=(0.0,findmax(PDOS1[50:end])[1]),
        # yticks = ([0],[L"$E_F$"]),
        )

    vline!(p, [0],
        lw=0.75,
        lc=LINE_COLOR[4],
        ls=:dashdot,
        label=false
        )
    # plot!(p, xlims=(fermi_energy-2, fermi_energy+2),ylims=(0, 7))
    annotate!(p, ANNOTATE_POSITION, 15, Plots.text("Fermi energy", ANNOTATE_SIZE, :dark, rotation = 90 ))
    return p
end

# plot_spin_pdos(FILE_DOS1,[["d"],["s","p"]],[[3,6,9,12],[51,54,57,60]],["up","down"])

# plot_spin_pdos(FILE_DOS1,[["d"],["s","p"]],[[13,16,19,22,25,28,31,34,37,40,43,46],collect(49:3:96)],["up","down"])