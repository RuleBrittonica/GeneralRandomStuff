using Plots
using Glob
using Base.Threads

# Parameters
L_min = 100.0    # Minimum rod length
L_max = 200.0    # Maximum rod length
L_step = 1    # Step size for rod lengths
Nx = 200         # Number of spatial points
Nt = 10000        # Number of time steps
Δt = 10          # Time step
α = 0.02         # Thermal diffusivity
T₀ = 14.76      # Initial temperature
Tc = 0.0        # Cold boundary condition
Th = 100.0      # Hot boundary condition

# Create an array of rod lengths
Ls = L_min:L_step:L_max

# Initialize the temperature matrix for each rod length
Ts = Dict{Float64, Matrix{Float64}}()

for L in Ls
    Δx = L / (Nx - 1)
    T = T₀ * ones(Nx, Nt+1)
    T[1, :] .= Th
    T[end, :] .= Tc
    Ts[L] = T
end

# Function to calculate the source term, which depends on position and time
function source_term(x, t, L)
    # Add your source term logic here. For demonstration, let's assume a simple case where the source
    # is a sinusoidal function of position and time.
    A = 10.0 # amplitude of the source term
    ω = 2π / L # frequency of the source term based on the length of the rod
    return A * sin(ω * x) * exp(-t)
end

# Create a threaded function to update temperature
function update_temperature_threaded(T, L, Δx, Δt, α, n)
    Threads.@threads for i in 2:Nx-1
        x_pos = (i-1) * Δx # position along the rod
        S = source_term(x_pos, n*Δt, L) # calculate the source term
        T[i, n+1] = T[i, n] + α * Δt / Δx^2 * (T[i+1, n] - 2*T[i, n] + T[i-1, n]) + S * Δt
    end
end

# Time-stepping loop for each rod
for (L, T) in Ts
    Δx = L / (Nx - 1)
    Threads.@threads for n in 1:Nt
        update_temperature_threaded(T, L, Δx, Δt, α, n)
    end
end

# Determine the maximum rod length for position array
max_length = maximum(Ls)
Δx_max = max_length / (Nx - 1)

# Create a vector of positions for the x-axis (common for all rods)
positions = collect(0:Δx_max:max_length)

# Initialize an empty array for the z values
z = Array{Float64}(undef, length(Ls), length(positions))

# Start the animation
heatmap_anim = @animate for n = 1:Nt
    # Construct each frame of the animation
    for (index, L) in enumerate(Ls)
        Δx = L / (Nx - 1)
        row = Ts[L][:, n]
        # Pad row to match the size of the longest rod if necessary
        if length(row) < length(positions)
            row = [row; zeros(length(positions) - length(row))]
        end
        z[index, :] = row
    end
    heatmap(positions, Ls, z, color=:thermal, clim=(0, Th), xlabel="Position (m)", ylabel="Rod Length (m)", title="Time: $(round(n*Δt, digits=0)) s")
end

# Get a list of existing GIFs in the output directory
existing_gifs = glob("thermal_rod_evolution_*.gif")

# Determine the next filename with an incremented number
next_number = isempty(existing_gifs) ? 1 : maximum([parse(Int, split(splitext(gif_file)[1], "_")[end]) for gif_file in existing_gifs]) + 1
output_filename = "thermal_rod_evolution_$(next_number).gif"

# Save the animation with the new filename
gif(heatmap_anim, output_filename, fps=15)

