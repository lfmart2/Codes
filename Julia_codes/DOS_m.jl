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

const FILE_FE = raw"/home/r_floren/BandStructure/QE/NbP/Main_folder/supercell4x4_2/NbP.supercell4x4.scf.out"
const FILE_DOS = raw"/home/r_floren/BandStructure/QE/NbP/github/Outputs/fire_cluster/NbP_spin/NbP_pdos.dat.pdos_tot"
const FILE_FOLDER = raw"/home/r_floren/BandStructure/QE/NbP/Main_folder/supercell4x4_2/"



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


    
function Band_structure_plot(FOLDER_PATH_BS::String,FILE_PATH_FE::String)
    p1 = plot();
    # Parameters for the plot
    LINE_WIDTH = 1.75
    LINE_COLOR = [:blue, :black]
    LINE_ALPHA = 0.5
    ANNOTATE_POSITION = 0.2
    ANNOTATE_SIZE = 8
    fermi_energy = extract_fermi_energy(FILE_PATH_FE)
    # Read bandstructure files from FOLDER_PATH_BS
    files_BS = filter(x->occursin(r"NbP\.supercell4x4_.\.bands\.dat",x), readdir(FOLDER_PATH_BS))
    x_axis = 0.0; x_label = zeros(length(files_BS)+1);
    for i in files_BS
        data = readdlm(FOLDER_PATH_BS*i)
        vc_size = findall(x->x==0.0, data[:,1])
        for j in 1:length(vc_size)-1
            plot!(p1,x_axis.+data[vc_size[j]:vc_size[j+1]-1,1],data[vc_size[j]:vc_size[j+1]-1,2],
                lw=LINE_WIDTH,
                lc=LINE_COLOR[1],
                ls=:dot,label=false,
                xticks = false)
        end
        plot!(p1,x_axis.+data[vc_size[end]:end-1,1],data[vc_size[end]:end-1,2],
            lw=LINE_WIDTH,
            lc=LINE_COLOR[2],
            ls=:dot,xticks = false,
            label=false)
        x_label[findall(x->x==i, files_BS)[1]+1] = data[vc_size[end]-1,1]
        x_axis += data[vc_size[end]-1,1]
    end
    hline!(p1,[fermi_energy],lw=0.75,lc=LINE_COLOR[2],ls=:dashdot,label=false)
    plot!(p1,xlims=(0,x_axis),
        grid=false,
        ylabel="Energy (eV)")
    vline!(p1,cumsum(x_label),
        lw=1,
        lc=:black,
        ls=:solid,
        la=LINE_ALPHA,
        label=false)
    plot!(p1,ylim=(4.0,5.2),
        xticks = (cumsum(x_label),["\$ Y \$","\$ \\Gamma \$","\$ X \$", "\$ M \$", "\$ \\Gamma \$"]))
    annotate!(p1,1.5,
        fermi_energy+ANNOTATE_POSITION,
        Plots.text("Fermi energy", ANNOTATE_SIZE))
    display(p1)
end

