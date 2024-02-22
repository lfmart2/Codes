using Plots
using DelimitedFiles
using CSV
using PeriodicTable

# # This code reads a file, scales the atomic positions by a
# # given factor, and writes the modified data back to the file.


const FILE_PWSCF = raw"/home/r_floren/BandStructure/QE/NbP/Local/supercell/NbP.supercell2x2.tst.in"
const ATOMIC_POSITIONS = "ATOMIC_POSITIONS"
const K_POINTS = "K_POINTS"
const SCALING_FACTOR = 0.8

function openfile(data::String)
    open(data, "r") do file
        return readlines(file)
    end
end

# Function to scale atomic positions in a file
function scale_atoms(FILE_PATH::String, scaling_factor::Float64)
    # Read the file into an array of lines
    lines = openfile(FILE_PATH)
    
    # Find the index of the line containing "ATOMIC_POSITIONS"
    atomic_positions_index = findfirst(i -> occursin(ATOMIC_POSITIONS, i), lines)
    
    # Find the index of the line containing "K_POINTS"
    k_points_index = findfirst(i -> occursin(K_POINTS, i), lines)

    # If either "ATOMIC_POSITIONS" or "K_POINTS" was not found, throw an error
    if atomic_positions_index isa Nothing || k_points_index isa Nothing
        error("Required keywords not found in the input data.")
    end

    # Extract the atomic positions from the lines between "ATOMIC_POSITIONS" and "K_POINTS"
    atomic_positions = map(x -> x[4], [split(lines[x]) for x in atomic_positions_index + 1:k_points_index - 1])
    
    # Scale the atomic positions by the scaling factor
    scaled_positions = parse.(Float64, atomic_positions) ./ scaling_factor
    
    # Split the lines containing the atomic positions into arrays of words
    positions_to_change = split.(lines[atomic_positions_index + 1:k_points_index - 1])

    # Replace the original atomic positions with the scaled positions
    for i in 1:length(positions_to_change)
        positions_to_change[i][4] = string(scaled_positions[i])
    end

    # Join the words back into lines and replace the original lines with the new lines
    lines[atomic_positions_index + 1:k_points_index - 1] .= "  " .* join.(positions_to_change, "  ")
    
    # Scale the last line of the file by the scaling factor
    lines[end] = "   " * join(string.(parse.(Float64, split(lines[end])) .* scaling_factor), "    ")

    # Return the modified lines
    return lines
end

# Call the function to scale the atomic positions in the file
scaled_lines = scale_atoms(FILE_PWSCF, SCALING_FACTOR)

# Write the modified lines back to the file
open(FILE_PWSCF, "w") do file
    println(file, join(scaled_lines, '\n'))
end