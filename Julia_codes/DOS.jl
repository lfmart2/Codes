using Plots
using DelimitedFiles
using CSV
using Plots.PlotMeasures



## Files to lead the QE output data
# const FILE_FE = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/fire_cluster/W.16.Bulk/nscf.out"
# const FILE_DOS = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/fire_cluster/NbP_spin/NbP_pdos.dat.pdos_tot"
# const FILE_SYM = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/Local/Bulk/NbP_band.labelinfo.dat"
# const FILE_BAND = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/Local/Bulk/bands.scf.dat.gnu"
# const FILE_BAND1 = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/fire_cluster/W.16.Bulk/NbP_band.dat"

const FILE_DOS1 = raw"/media/r_floren/FernandoM/Bandstructure/VASP/NbP/Supercell_final/Coord_1/nSpin_H/DOS/DOSCAR"
const FILE_DOS2 = raw"/home/r_floren/BandStructure/VASP/NbP/Supercell_2_2_3/non-spin-QE/DOS2/DOSCAR"
const FILE_DOS5 = raw"/home/r_floren/BandStructure/VASP/NbP/Supercell_final/Coord_1/H/out/DOSCAR"




function openfile(data::String)
    open(data, "r") do file
        return readlines(file)
    end
end

function extract_fermi_energy(FILE_PATH::String)::Float64
    return parse(Float64, split(openfile(FILE_PATH)[6])[4])
end


function plot_dos(FILE_PATH_DOSCAR::String,e_min::Float64,e_max::Float64)
    # Range of Energy
    #E_min = readdlm(FILE_PATH_DOSCAR)[6,2];
    #E_max = readdlm(FILE_PATH_DOSCAR)[6,1];
    # Parameters for the plot
    LINE_WIDTH = 1.75
    LINE_COLOR = [:black, :black]
    FILL_ALPHA = 0.2
    ANNOTATE_POSITION = -0.05
    ANNOTATE_SIZE = 12
    fermi_energy = extract_fermi_energy(FILE_PATH_DOSCAR)
    e_mesh = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[6])[3])
    energy_data = readdlm(FILE_PATH_DOSCAR)[7:e_mesh+6, 1:3]
    p = plot(leftmargin = 5mm)
    tmp = findmin(abs.(energy_data[1:end,1] .- fermi_energy))[2]
    plot!(p, energy_data[1:end,1].-fermi_energy, energy_data[1:end,2],
        lw=LINE_WIDTH,
        lc=LINE_COLOR[1],
        label=false
        )
    plot!(p, energy_data[1:tmp,1].-fermi_energy, energy_data[1:tmp,2],
        fill = (FILL_ALPHA, 0, LINE_COLOR[1]),
        lw=LINE_WIDTH,
        lc=LINE_COLOR[1],
        label=false
        )
    plot!(p, grid=false, yticks=false,
        xlabel="Energy (eV)",
        ylabel="DOS (states/eV)",
        xlims=(e_min, e_max),
        ylims=(0.0,findmax(energy_data[1:end,2])[1]+100)
        )
    plot!(p, fermi_energy*ones(2), [0, findmax(energy_data[1:end,2])[1]],
        lw=0.75,
        lc=LINE_COLOR[2],
        ls=:dashdot,
        label=false
        )
    # plot!(p, xlims=(fermi_energy-2, fermi_energy+2),ylims=(0, 7))
    annotate!(p, fermi_energy + ANNOTATE_POSITION, 150, Plots.text("Fermi energy", ANNOTATE_SIZE, :dark, rotation = 90 ))
    return p
end

function plot_spin_dos(FILE_PATH_DOSCAR::String,e_min::Float64,e_max::Float64)
    # Range of Energy
    #E_min = readdlm(FILE_PATH_DOSCAR)[6,2];
    #E_max = readdlm(FILE_PATH_DOSCAR)[6,1];
    # Parameters for the plot
    LINE_WIDTH = 1.75
    LINE_COLOR = [:red, :blue, :gray]
    FILL_ALPHA = 0.2
    ANNOTATE_POSITION = -0.4
    ANNOTATE_SIZE = 12
    fermi_energy = extract_fermi_energy(FILE_PATH_DOSCAR)
    e_mesh = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[6])[3])
    energy_data = readdlm(FILE_PATH_DOSCAR)[7:e_mesh+6, 1:3]
    p = plot(leftmargin = 5mm)
    plot!(p, energy_data[1:end,2], energy_data[1:end,1].-fermi_energy,
        lw=LINE_WIDTH,
        lc=LINE_COLOR[1],
        label  = false
        # label="spin-up"
        )
    plot!(p, energy_data[1:end,3], energy_data[1:end,1].-fermi_energy,
        lw=LINE_WIDTH,
        lc=LINE_COLOR[2],
        label = false
        # label="spin-down"
        )
    energy_data[findmin(abs.(energy_data[1:end,1] .- fermi_energy))[2],2] = 0.0;
    energy_data[findmin(abs.(energy_data[1:end,1] .- fermi_energy))[2],3] = 0.0;
    plot!(p, energy_data[1:findmin(abs.(energy_data[1:end,1] .- fermi_energy))[2],2], energy_data[1:findmin(abs.(energy_data[1:end,1] .- fermi_energy))[2],1].-fermi_energy,
        fillalpha = FILL_ALPHA,
        fillrange = -20,
        fillcolor = LINE_COLOR[3],
        lw        = 0,
        lc        = LINE_COLOR[1],
        label     = false
        )
    plot!(p, energy_data[1:findmin(abs.(energy_data[1:end,1] .- fermi_energy))[2],3], energy_data[1:findmin(abs.(energy_data[1:end,1] .- fermi_energy))[2],1].-fermi_energy,
        fillalpha = FILL_ALPHA,
        fillrange = -20,
        fillcolor = LINE_COLOR[3],
        lw        = 0,
        lc        = LINE_COLOR[2],
        label     = false,
        )
    plot!(p, grid=false, xticks=false, yticks=false,
        # ylabel="Energy (eV)",
        # xlabel="DOS (states/eV)",
        ylims=(fermi_energy-e_min, fermi_energy+e_max),
        xlims=(0.0,findmax(energy_data[10:end,2])[1])
        )
    plot!(p, [findmax(energy_data[1:end,2])[1], 0], 0.0*ones(2),
        lw=0.75,
        lc=LINE_COLOR[3],
        ls=:dashdot,
        label=false
        )
    # plot!(p, xlims=(fermi_energy-2, fermi_energy+2),ylims=(0, 7))
    # annotate!(p, fermi_energy + ANNOTATE_POSITION, 150, Plots.text("Fermi energy", ANNOTATE_SIZE, :dark, rotation = 90 ))
    return p
end

function plot_spin_dos2(FILE_PATH_DOSCAR::String,e_min::Float64,e_max::Float64)
    # Range of Energy
    #E_min = readdlm(FILE_PATH_DOSCAR)[6,2];
    #E_max = readdlm(FILE_PATH_DOSCAR)[6,1];
    # Parameters for the plot
    LINE_WIDTH = 1.5
    LINE_COLOR = [:red, :blue, :black]
    FILL_ALPHA = 0.5
    ANNOTATE_POSITION = -0.05
    ANNOTATE_SIZE = 12
    fermi_energy = extract_fermi_energy(FILE_PATH_DOSCAR)
    e_mesh = parse(Int64, split(openfile(FILE_PATH_DOSCAR)[6])[3])
    energy_data = readdlm(FILE_PATH_DOSCAR)[7:e_mesh+6, 1:3]
    p = plot(leftmargin = 5mm)
    tmp = findmin(abs.(energy_data[1:end,1] .- fermi_energy))[2]
    plot!(p, energy_data[1:tmp,1],energy_data[1:tmp,2].+ energy_data[1:tmp,3],
        fill = (FILL_ALPHA, 0, LINE_COLOR[3]),
        lw=LINE_WIDTH,
        lc=LINE_COLOR[3],
        label=false
    )
    plot!(p, energy_data[1:end,1], energy_data[1:end,2] .+ energy_data[1:end,3],
        lw=LINE_WIDTH,
        lc=LINE_COLOR[3],
        label="Total"
        )
    plot!(p, energy_data[1:end,1], energy_data[1:end,2] ,
        lw=LINE_WIDTH,
        lc=LINE_COLOR[1],
        label="spin-up"
        )
    plot!(p, energy_data[1:end,1], energy_data[1:end,3],
        lw=LINE_WIDTH,
        lc=LINE_COLOR[2],
        label="spin-down"
        )
    plot!(p, grid=false, yticks=false,
        xlabel="Energy (eV)",
        ylabel="DOS (states/eV)",
        xlims=(fermi_energy-e_min,fermi_energy+ e_max),
        ylims=(0.0,findmax(energy_data[60:end,2])[1])
        )
    plot!(p, fermi_energy*ones(2), [0, findmax(energy_data[1:end,2])[1]],
        lw=0.75,
        lc=LINE_COLOR[3],
        ls=:dashdot,
        label=false
        )
    # plot!(p, xlims=(fermi_energy-2, fermi_energy+2),ylims=(0, 7))
    annotate!(p, fermi_energy + ANNOTATE_POSITION, 150, Plots.text("Fermi energy", ANNOTATE_SIZE, :dark, rotation = 90 ))
    return p
end