# include the jfdtd package from the direcotry
include("src/jfdtd.jl")
using .jfdtd

fdtd_2d_guas_pulse_nopml(400)

fdtd_2d_guas_pulse_pml(20, 8)

fdtd_2d_planewave_pml(200, 8)


# Run the simulation with sample inputs
nsteps = 400
radius = 3.0
epsilon = 4.0
sigma = 0.0

p = fdtd_2d_ep_sphere_nopml(nsteps, radius, epsilon, sigma)

# Run the simulation with the option to save the GIF
nsteps = 200  # Adjust as needed
fdtd_3d_dipole_nopml(nsteps, save_gif = true)

# Run the simulation
nsteps = 240  # Adjust as needed
fdtd_3d_sphere_pml(nsteps, save_gif = true)

