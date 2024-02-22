using Plots
using DelimitedFiles
using CSV

const LINE_COLOR = [:red, :green, :blue]
const LINE_WIDTH = 2
const LINE_ALPHA = 0.5

function plot_line(result, i, j, color)
    plot!([result[i,1],result[i,1]+result[j,1]], [result[i,2],result[i,2]+result[j,2]],
    [result[i,3],result[i,3]+result[j,3]], labels = false, lc = color, lw=LINE_WIDTH, la=LINE_ALPHA)
end

function plot_unit_cell(data_file::String)
    datafile = open(data_file,"r")
    lines = readlines(datafile)
    match = findfirst(i -> occursin("CELL_PARAMETERS", i), lines)
    cell_parameters = match isa Nothing ? nothing : hcat([parse.(Float64, split(lines[i])) for i in match+1:match+3]...)'

    for i in 1:size(cell_parameters)[1]
        plot_line(cell_parameters, i, i, LINE_COLOR[i])
    end

    for i in 1:3
        for j in i+1:3
            plot_line(cell_parameters, i, j, :black)
            plot_line(cell_parameters, j, i, :black)
        end
    end

    plot_line(cell_parameters, 1, 2, :black)
    plot_line(cell_parameters, 2, 3, :black)
    plot_line(cell_parameters, 1, 3, :black)
end