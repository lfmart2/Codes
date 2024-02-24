using Plots
using DelimitedFiles
using CSV
using Plots.PlotMeasures


const LINE_WIDTH = 0.75
const LINE_COLOR = [:red, :black]
const FILL_ALPHA = 0.5
const ANNOTATE_POSITION = -0.2
const ANNOTATE_SIZE = 12

# Files to lead the QE output data
const FILE_FE = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/Local/Bulk/nscf.out"
const FILE_DOS = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/fire_cluster/NbP_spin/NbP_pdos.dat.pdos_tot"
const FILE_SYM = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/Local/Bulk/NbP_band.labelinfo.dat"
const FILE_BAND = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/Local/Bulk/bands.scf.dat.gnu"
const FILE_BAND1 = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/Local/Bulk/NbP_band.dat"


function openfile(data::String)
    open(data, "r") do file
        return readlines(file)
    end
end

function extract_fermi_energy(FILE_PATH::String)::Float64
    lines = openfile(FILE_PATH)
    match = findfirst(i -> occursin("Fermi", i), lines)
    if match isa Nothing
        error("Fermi energy not found in the input data.")
    else
        return parse(Float64, split(lines[match])[end-1])
    end
end

const FERMI_ENERGY = extract_fermi_energy(FILE_FE)

function plot_dos(FILE_PATH_DOS::String, FILE_PATH_FE::String)
    fermi_energy = extract_fermi_energy(FILE_PATH_FE)
    energy_data = readdlm(FILE_PATH_DOS)
    p = plot(leftmargin = 5mm)
    plot!(p, energy_data[2:end,1], energy_data[2:end,2], lw=LINE_WIDTH, lc=LINE_COLOR[1], label=false)
    plot!(p, energy_data[2:findmin(abs.(energy_data[2:end,1] .- fermi_energy))[2],1], energy_data[2:findmin(abs.(energy_data[2:end,1] .- fermi_energy))[2],2], fill = (FILL_ALPHA, 0, LINE_COLOR[1]), lw=LINE_WIDTH, lc=LINE_COLOR[1], label=false)
    plot!(p, grid=false, yticks=false, xlabel="Energy (eV)", ylabel="DOS (states/eV)", xlims=(energy_data[2,1], energy_data[end,1]), ylims=(energy_data[2,2], findmax(energy_data[2:end,2])[1]))
    plot!(p, fermi_energy*ones(2), [0, findmax(energy_data[2:end,2])[1]], lw=LINE_WIDTH, lc=LINE_COLOR[2], ls=:dashdot, label=false)
    plot!(p, xlims=(8, 19),ylims=(0, 7))
    annotate!(p, fermi_energy + ANNOTATE_POSITION, 5, Plots.text("Fermi energy", ANNOTATE_SIZE, :dark, rotation = 90 ))
    display(p)
end

function sym_points(FILE_PATH::String)
    sym_points = openfile(FILE_PATH)
    indices = findall(i -> occursin("high-symmetry point:", i), sym_points)
    vctr = [parse(Float64, split(sym_points[i])[end]) for i in indices]
    return unique(vctr)
end


    
function Band_structure_plot(FILE_PATH_BS::String,FILE_PATH_SYM::String,FILE_PATH_FE::String)
    p1 = plot();
    data = readdlm(FILE_PATH_BS)
    fermi_energy = extract_fermi_energy(FILE_PATH_FE)
    vc_size = findall(x->x==0.0, data[:,1])
    for i in 1:length(vc_size)-1
        plot!(p1,data[vc_size[i]:vc_size[i+1]-1,1]./data[end,1],data[vc_size[i]:vc_size[i+1]-1,2],lw=0.75,lc=:red,ls=:dash,label=false,xticks = false)
    end
    plot!(p1,data[vc_size[end]:end-1,1]./data[end,1],data[vc_size[end]:end-1,2],lw=0.75,lc=:red,ls=:dash,label="SCF",xticks = false)
    hline!(p1,[fermi_energy],lw=0.75,lc=:black,ls=:dashdot,label=false)
    vline!(p1,readdlm(FILE_PATH_SYM)[:,3]./readdlm(FILE_PATH_SYM)[end,3], lw=1,lc=:black,ls=:solid,label=false)
    plot!(p1,xlims=(0,1),grid=false, ylabel="Energy (eV)")
    plot!(p1,ylim=(8,18.5), xticks = (readdlm(FILE_PATH_SYM)[:,3]./readdlm(FILE_PATH_SYM)[end,3],["\$\\Gamma \$","\$X \$","\$ M \$","\$ Y \$","\$ K \$","\$ F \$","\$ \\Gamma \$"]))
    annotate!(p1,1.5,fermi_energy+0.2, Plots.text("Fermi energy", 8))
    display(p1)
end

function Band_structure_plot1(FILE_PATH_BS::String,FILE_PATH_SYM::String,FILE_PATH_FE::String)
    data = readdlm(FILE_PATH_BS)
    vc_size = findall(x->x==0.0, data[:,1])
    for i in 1:length(vc_size)-1
        plot!(data[vc_size[i]:vc_size[i+1]-1,1]./data[end,1],data[vc_size[i]:vc_size[i+1]-1,2],lw=0.75,lc=:blue,label=false)
    end
    plot!(data[vc_size[end]:end-1,1]./data[end,1],data[vc_size[end]:end-1,2],lw=0.75,lc=:blue,label="Tight-Binding")
end
