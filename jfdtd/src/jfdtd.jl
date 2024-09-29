module jfdtd

"""
Working on 1D, 2D, and 3D FDTD algorithms using julia. 
Initially based on Sullivan's book on FDTD:
[1] D. M. Sullivan, Electromagnetic Simulation Using The FDTD Method. 
    IEEE Press Series on Electromagnetic Wave Theory. Wiley-IEEE Press, 2013.
"""

using Plots
export fdtd_2d_guas_pulse_nopml, fdtd_2d_guas_pulse_pml, fdtd_2d_planewave_pml, fdtd_2d_ep_sphere_nopml, fdtd_3d_dipole_nopml, fdtd_3d_sphere_pml
####################################################


function fdtd_2d_guas_pulse_nopml(nsteps_input)
    """
    Basic 2D TM FDTD Simulation with Animation in Julia no plots guassian pulse source

    """
    IE = 60
    JE = 60

    # Initialize arrays
    ga = ones(Float64, IE, JE)
    dz = zeros(Float64, IE, JE)
    ez = zeros(Float64, IE, JE)
    hx = zeros(Float64, IE, JE)
    hy = zeros(Float64, IE, JE)

    # Simulation parameters
    ic = div(IE, 2)
    jc = div(JE, 2)
    ddx = 0.01          # Cell size (meters)
    dt = ddx / (6e8)    # Time step (seconds)
    epsz = 8.8e-12      # Permittivity of free space
    pi = 3.14159

    # Source parameters
    t0 = 20.0
    spread = 6.0
    T = 0.0

    nsteps = nsteps_input

    # Initialize animation
    anim = @animate for n = 1:nsteps
        T += 1.0

        # Update Dz field
        for j = 2:JE-1
            for i = 2:IE-1
                dz[i,j] += 0.5 * (hy[i,j] - hy[i-1,j] - hx[i,j] + hx[i,j-1])
            end
        end

        # Gaussian pulse source in the middle
        pulse = exp(-0.5 * ((t0 - T) / spread)^2)
        dz[ic,jc] = pulse

        # Update Ez field
        for j = 2:JE-1
            for i = 2:IE-1
                ez[i,j] = ga[i,j] * dz[i,j]
            end
        end

        # Update Hx field
        for j = 1:JE-1
            for i = 1:IE-1
                hx[i,j] += 0.5 * (ez[i,j] - ez[i,j+1])
            end
        end

        # Update Hy field
        for j = 1:JE-1
            for i = 1:IE-1
                hy[i,j] += 0.5 * (ez[i+1,j] - ez[i,j])
            end
        end

        # Visualization
        heatmap(
            ez',
            title = "Ez Field at T = $T",
            clims = (-1, 1),
            xlabel = "X",
            ylabel = "Y",
            colorbar = false,
            framestyle = :none
        )
    end

    # Save animation as GIF
    gif(anim, "fd2d_3_1_animation.gif", fps = 20)
end



function fdtd_2d_guas_pulse_pml(nsteps_input, npml_input)
    """
    2D TM FDTD Simulation with PML in Julia with guassian pulse source

    """
    IE = 60
    JE = 60

    # Initialize arrays
    ga = ones(Float64, IE, JE)
    dz = zeros(Float64, IE, JE)
    ez = zeros(Float64, IE, JE)
    hx = zeros(Float64, IE, JE)
    hy = zeros(Float64, IE, JE)
    ihx = zeros(Float64, IE, JE)
    ihy = zeros(Float64, IE, JE)

    # Simulation parameters
    ic = div(IE, 2) - 10
    jc = div(JE, 2) - 10
    ddx = 0.01          # Cell size
    dt = ddx / (6e8)    # Time step
    epsz = 8.8e-12
    pi = 3.14159

    # PML parameters
    npml = npml_input
    gi2 = ones(Float64, IE)
    gi3 = ones(Float64, IE)
    fi1 = zeros(Float64, IE)
    fi2 = ones(Float64, IE)
    fi3 = ones(Float64, IE)
    gj2 = ones(Float64, JE)
    gj3 = ones(Float64, JE)
    fj1 = zeros(Float64, JE)
    fj2 = ones(Float64, JE)
    fj3 = ones(Float64, JE)

    # Calculate PML parameters
    for i in 1:npml
        xnum = npml - i + 0.5
        xd = npml
        xxn = xnum / xd
        xn = 0.25 * xxn^3

        gi2[i] = 1.0 / (1.0 + xn)
        gi2[IE - i + 1] = gi2[i]
        gi3[i] = (1.0 - xn) / (1.0 + xn)
        gi3[IE - i + 1] = gi3[i]

        fi1[i] = xn
        fi1[IE - i] = xn
        fi2[i] = 1.0 / (1.0 + xn)
        fi2[IE - i] = fi2[i]
        fi3[i] = (1.0 - xn) / (1.0 + xn)
        fi3[IE - i] = fi3[i]
    end

    for j in 1:npml
        xnum = npml - j + 0.5
        xd = npml
        xxn = xnum / xd
        xn = 0.25 * xxn^3

        gj2[j] = 1.0 / (1.0 + xn)
        gj2[JE - j + 1] = gj2[j]
        gj3[j] = (1.0 - xn) / (1.0 + xn)
        gj3[JE - j + 1] = gj3[j]

        fj1[j] = xn
        fj1[JE - j] = xn
        fj2[j] = 1.0 / (1.0 + xn)
        fj2[JE - j] = fj2[j]
        fj3[j] = (1.0 - xn) / (1.0 + xn)
        fj3[JE - j] = fj3[j]
    end

    # Source parameters
    t0 = 40.0
    spread = 12.0
    T = 0.0

    nsteps = nsteps_input

    # Initialize animation
    anim = @animate for n = 1:nsteps
        T += 1.0

        # Update Dz field
        for j = 2:JE-1
            for i = 2:IE-1
                dz[i,j] = gi3[i]*gj3[j]*dz[i,j] + gi2[i]*gj2[j]*0.5*(hy[i,j] - hy[i-1,j] - hx[i,j] + hx[i,j-1])
            end
        end

        # Gaussian pulse source
        pulse = exp(-0.5 * ((T - t0) / spread)^2)
        dz[ic,jc] = pulse

        # Update Ez field
        for j = 2:JE-1
            for i = 2:IE-1
                ez[i,j] = ga[i,j]*dz[i,j]
            end
        end

        # Update Hx field
        for j = 1:JE-1
            for i = 1:IE
                curl_e = ez[i,j] - ez[i,j+1]
                ihx[i,j] += fi1[i]*curl_e
                hx[i,j] = fj3[j]*hx[i,j] + fj2[j]*0.5*(curl_e + ihx[i,j])
            end
        end

        # Update Hy field
        for j = 1:JE
            for i = 1:IE-1
                curl_e = ez[i+1,j] - ez[i,j]
                ihy[i,j] += fj1[j]*curl_e
                hy[i,j] = fi3[i]*hy[i,j] + fi2[i]*0.5*(curl_e + ihy[i,j])
            end
        end

        # Visualization
        heatmap(
            ez',
            title = "Ez Field with PML at T = $T",
            clims = (-1, 1),
            xlabel = "X",
            ylabel = "Y",
            colorbar = false,
            framestyle = :none
        )
    end

    # Save animation as GIF
    gif(anim, "fd2d_3_2_animation.gif", fps = 20)
end

function fdtd_2d_planewave_pml(nsteps_input, npml_input)
    """
    2D TM FDTD Simulation with PML in Julia with plane wave source

    """
    IE = 60
    JE = 60

    # Initialize arrays
    ga = ones(Float64, IE, JE)
    dz = zeros(Float64, IE, JE)
    ez = zeros(Float64, IE, JE)
    hx = zeros(Float64, IE, JE)
    hy = zeros(Float64, IE, JE)
    ihx = zeros(Float64, IE, JE)
    ihy = zeros(Float64, IE, JE)
    ez_inc = zeros(Float64, JE+2)  # Increased size by 1
    hx_inc = zeros(Float64, JE+1)

    ez_inc_low_m1 = 0.0
    ez_inc_low_m2 = 0.0
    ez_inc_high_m1 = 0.0
    ez_inc_high_m2 = 0.0

    # Simulation parameters
    ic = div(IE, 2)
    jc = div(JE, 2)

    ia = 7  # Total/scattered field boundaries
    ib = IE - ia - 1
    ja = 7
    jb = JE - ja - 1

    ddx = 0.01          # Cell size
    dt = ddx / (6e8)    # Time step
    epsz = 8.8e-12
    pi = 3.14159

    # PML parameters
    npml = npml_input
    gi2 = ones(Float64, IE)
    gi3 = ones(Float64, IE)
    fi1 = zeros(Float64, IE)
    fi2 = ones(Float64, IE)
    fi3 = ones(Float64, IE)
    gj2 = ones(Float64, JE)
    gj3 = ones(Float64, JE)
    fj1 = zeros(Float64, JE)
    fj2 = ones(Float64, JE)
    fj3 = ones(Float64, JE)

    # Calculate PML parameters
    for i in 1:npml
        xnum = npml - i + 0.5
        xd = npml
        xxn = xnum / xd
        xn = 0.25 * xxn^3

        gi2[i] = 1.0 / (1.0 + xn)
        gi2[IE - i + 1] = gi2[i]
        gi3[i] = (1.0 - xn) / (1.0 + xn)
        gi3[IE - i + 1] = gi3[i]

        fi1[i] = xn
        fi1[IE - i] = xn
        fi2[i] = 1.0 / (1.0 + xn)
        fi2[IE - i] = fi2[i]
        fi3[i] = (1.0 - xn) / (1.0 + xn)
        fi3[IE - i] = fi3[i]
    end

    for j in 1:npml
        xnum = npml - j + 0.5
        xd = npml
        xxn = xnum / xd
        xn = 0.25 * xxn^3

        gj2[j] = 1.0 / (1.0 + xn)
        gj2[JE - j + 1] = gj2[j]
        gj3[j] = (1.0 - xn) / (1.0 + xn)
        gj3[JE - j + 1] = gj3[j]

        fj1[j] = xn
        fj1[JE - j] = xn
        fj2[j] = 1.0 / (1.0 + xn)
        fj2[JE - j] = fj2[j]
        fj3[j] = (1.0 - xn) / (1.0 + xn)
        fj3[JE - j] = fj3[j]
    end

    # Source parameters
    t0 = 20.0
    spread = 8.0
    T = 0.0

    nsteps = nsteps_input

    # Initialize animation
    anim = @animate for n = 1:nsteps
        T += 1.0

        # Update incident Ez field
        for j = 2:JE
            ez_inc[j] += 0.5 * (hx_inc[j-1] - hx_inc[j])
        end

        # ABC for incident field
        ez_inc[1] = ez_inc_low_m2
        ez_inc_low_m2 = ez_inc_low_m1
        ez_inc_low_m1 = ez_inc[2]

        ez_inc[JE+1] = ez_inc_high_m2
        ez_inc_high_m2 = ez_inc_high_m1
        ez_inc_high_m1 = ez_inc[JE]

        # Update Dz field
        for j = 2:JE-1
            for i = 2:IE-1
                dz[i,j] = gi3[i]*gj3[j]*dz[i,j] + gi2[i]*gj2[j]*0.5*(hy[i,j] - hy[i-1,j] - hx[i,j] + hx[i,j-1])
            end
        end

        # Plane wave source
        pulse = exp(-0.5 * ((t0 - T) / spread)^2)
        ez_inc[4] = pulse

        # Incident Dz values at the boundaries
        for i = ia:ib
            dz[i,ja] += 0.5 * hx_inc[ja-1]
            dz[i,jb] -= 0.5 * hx_inc[jb]
        end

        # Update Ez field
        for j = 2:JE-1
            for i = 2:IE-1
                ez[i,j] = ga[i,j]*dz[i,j]
            end
        end

        # Update incident Hx field
        for j = 1:JE
            hx_inc[j] += 0.5 * (ez_inc[j] - ez_inc[j+1])
        end

        # Update Hx field
        for j = 1:JE-1
            for i = 1:IE
                curl_e = ez[i,j] - ez[i,j+1]
                ihx[i,j] += fi1[i]*curl_e
                hx[i,j] = fj3[j]*hx[i,j] + fj2[j]*0.5*(curl_e + ihx[i,j])
            end
        end

        # Incident Hx values at the boundaries
        for i = ia:ib
            hx[i,ja-1] += 0.5 * ez_inc[ja]
            hx[i,jb]   -= 0.5 * ez_inc[jb]
        end

        # Update Hy field
        for j = 1:JE
            for i = 1:IE-1
                curl_e = ez[i+1,j] - ez[i,j]
                ihy[i,j] += fj1[j]*curl_e
                hy[i,j] = fi3[i]*hy[i,j] + fi2[i]*0.5*(curl_e + ihy[i,j])
            end
        end

        # Incident Hy values at the boundaries
        for j = ja:jb
            hy[ia-1,j] -= 0.5 * ez_inc[j]
            hy[ib,j]   += 0.5 * ez_inc[j]
        end

        # Visualization
        heatmap(
            ez',
            title = "Ez Field with Plane Wave at T = $T",
            clims = (0, 1),
            xlabel = "X",
            ylabel = "Y",
            colorbar = true,
            framestyle = :none,
            color = :RdBu,
        )
    end

    # Save animation as GIF
    gif(anim, "fd2d_3_3_animation.gif", fps = 20)
end


function fdtd_2d_ep_sphere_nopml(nsteps_input::Int64, radius_input::Float64, epsilon_input::Float64, sigma_input::Float64)
    """
    2D TM FDTD Simulation with a dielectric sphere in Julia with no PML

    """
    # Physical constants
    epsz = 8.854e-12              # Vacuum permittivity (F/m)
    μ0 = 4π * 1e-7                # Vacuum permeability (H/m)
    c0 = 1 / sqrt(epsz * μ0)      # Speed of light in vacuum (m/s)
    pi = π                        # Pi constant

    # Simulation parameters
    IE = 100                       # Number of grid points in x-direction
    JE = 100                       # Number of grid points in y-direction
    dx = 0.005                     # Spatial step size (m)
    dy = dx                       # For simplicity, dx = dy
    dt = dx / (2 * c0)            # Time step size (s), satisfying the Courant condition
    nsteps = nsteps_input         # Number of time steps
    radius = radius_input         # Radius of the cylinder (in cells)
    epsilon = epsilon_input       # Relative permittivity of the cylinder
    sigma = sigma_input           # Conductivity of the cylinder

    # Frequencies for Fourier Transform analysis
    NFREQS = 3
    freq = [50e6, 300e6, 700e6]           # Frequencies in Hz
    arg = 2 * pi * freq * dt              # Angular frequency times time step

    # Total field/scattered field boundaries
    ia = 7
    ib = IE - ia - 1
    ja = 7
    jb = JE - ja - 1

    # Initialize fields
    Dz = zeros(IE, JE)
    Ez = zeros(IE, JE)
    Ez_old = zeros(IE, JE)
    Hx = zeros(IE, JE)
    Hy = zeros(IE, JE)
    iz = zeros(IE, JE)
    ihx = zeros(IE, JE)
    ihy = zeros(IE, JE)

    # Incident fields for plane wave
    ez_inc = zeros(JE + 1)
    hx_inc = zeros(JE + 1)  # Adjusted size to prevent out-of-bounds access

    # Fourier Transform arrays
    real_pt = zeros(NFREQS, IE, JE)
    imag_pt = zeros(NFREQS, IE, JE)
    amp = zeros(IE, JE)
    phase = zeros(IE, JE)
    real_in = zeros(NFREQS)
    imag_in = zeros(NFREQS)
    amp_in = zeros(NFREQS)
    phase_in = zeros(NFREQS)

    # Material properties arrays
    ga = ones(IE, JE)
    gb = zeros(IE, JE)

    # Coefficient arrays
    gi2 = ones(IE)
    gi3 = ones(IE)
    fi1 = zeros(IE)
    fi2 = ones(IE)
    fi3 = ones(IE)
    gj2 = ones(JE)
    gj3 = ones(JE)
    fj1 = zeros(JE)
    fj2 = ones(JE)
    fj3 = ones(JE)

    # Define the dielectric cylinder
    ic = IE ÷ 2
    jc = JE ÷ 2

    for j in ja:jb
        for i in ia:ib
            xdist = (ic - i)
            ydist = (jc - j)
            dist = sqrt(xdist^2 + ydist^2)
            if dist <= radius
                ga[i, j] = 1.0 / (epsilon + (sigma * dt / epsz))
                gb[i, j] = sigma * dt / epsz
            end
        end
    end

    # Simulation constants
    t0 = 20.0
    spread = 6.0

    # Prepare the plot
    p = heatmap(
        Ez',
        clim = (-0.01, 0.01),
        aspect_ratio = 1,
        title = "Ez Field at Time Step 0",
        xlabel = "X",
        ylabel = "Y",
        color = :RdBu,
        xlims = (1, IE),
        ylims = (1, JE),
        framestyle = :box
    )

    # Time-stepping loop
    for n = 1:nsteps
        T = n

        # --- Update incident Ez field ---
        for j = 2:JE
            ez_inc[j] += 0.5 * (hx_inc[j - 1] - hx_inc[j])
        end

        # Fourier Transform of the incident field
        for m = 1:NFREQS
            real_in[m] += cos(arg[m] * T) * ez_inc[ja]
            imag_in[m] -= sin(arg[m] * T) * ez_inc[ja]
        end

        # --- Update Dz field ---
        for j = 2:JE - 1
            for i = 2:IE - 1
                Dz[i, j] = gi3[i] * gj3[j] * Dz[i, j] + gi2[i] * gj2[j] * 0.5 * (
                    Hy[i, j] - Hy[i - 1, j] - Hx[i, j] + Hx[i, j - 1]
                )
            end
        end

        # --- Apply source ---
        pulse = exp(-0.5 * ((t0 - T) / spread)^2)
        ez_inc[4] = pulse

        # --- Incident Dz values ---
        for i = ia:ib
            Dz[i, ja] += 0.5 * hx_inc[ja - 1]
            Dz[i, jb] -= 0.5 * hx_inc[jb]
        end

        # --- Update Ez field ---
        for j = 2:JE - 1
            for i = 2:IE - 1
                Ez[i, j] = ga[i, j] * (Dz[i, j] - iz[i, j])
                iz[i, j] += gb[i, j] * Ez[i, j]
            end
        end

        # --- Fourier Transform of Ez field ---
        for j = 1:JE
            for i = 1:IE
                for m = 1:NFREQS
                    real_pt[m, i, j] += cos(arg[m] * T) * Ez[i, j]
                    imag_pt[m, i, j] += sin(arg[m] * T) * Ez[i, j]
                end
            end
        end

        # --- Update incident Hx field ---
        for j = 1:JE
            hx_inc[j] += 0.5 * (ez_inc[j] - ez_inc[j + 1])
        end

        # --- Update Hx field ---
        for j = 1:JE - 1
            for i = 1:IE
                curl_e = Ez[i, j] - Ez[i, j + 1]
                ihx[i, j] += fi1[i] * curl_e
                Hx[i, j] = fj3[j] * Hx[i, j] + fj2[j] * 0.5 * (curl_e + ihx[i, j])
            end
        end

        # --- Incident Hx values ---
        for i = ia:ib
            Hx[i, ja - 1] += 0.5 * ez_inc[ja]
            Hx[i, jb] -= 0.5 * ez_inc[jb]
        end

        # --- Update Hy field ---
        for j = 1:JE
            for i = 1:IE - 1
                curl_e = Ez[i + 1, j] - Ez[i, j]
                ihy[i, j] += fj1[j] * curl_e
                Hy[i, j] = fi3[i] * Hy[i, j] + fi2[i] * 0.5 * (curl_e + ihy[i, j])
            end
        end

        # --- Incident Hy values ---
        for j = ja:jb
            Hy[ia - 1, j] -= 0.5 * ez_inc[j]
            Hy[ib, j] += 0.5 * ez_inc[j]
        end

        # --- Update the plot ---
        if n % 5 == 0 || n == 1
            # Update the data in the heatmap
            p.series_list[1][:z] = Ez'
            p.series_list[1][:title] = "Ez Field at Time Step $n"
            display(p)
            sleep(0.01)  # Pause briefly to allow the plot to update
        end
    end

    # Return the plot object
    return p
end


function fdtd_3d_dipole_nopml(nsteps_input::Int64; save_gif::Bool = false)
    """
    3D FDTD Simulation with a dipole source in Julia with no PML

    """
    # Physical constants
    epsz = 8.854e-12              # Vacuum permittivity (F/m)
    μ0 = 4π * 1e-7                # Vacuum permeability (H/m)
    c0 = 1 / sqrt(epsz * μ0)      # Speed of light in vacuum (m/s)

    # Simulation parameters
    IE = 100                       # Number of grid points in x-direction
    JE = 100                       # Number of grid points in y-direction
    KE = 100                       # Number of grid points in z-direction
    ddx = 0.01                    # Spatial step size (m)
    dt = ddx / (6e8)              # Time step size (s)
    nsteps = nsteps_input         # Number of time steps

    # Pulse parameters matching the C code
    t0 = 20.0
    spread = 6.0

    # Center indices
    ic = div(IE, 2)
    jc = div(JE, 2)
    kc = div(KE, 2)

    # Initialize fields
    ex = zeros(Float64, IE, JE, KE)
    ey = zeros(Float64, IE, JE, KE)
    ez = zeros(Float64, IE, JE, KE)
    dx = zeros(Float64, IE, JE, KE)
    dy = zeros(Float64, IE, JE, KE)
    dz = zeros(Float64, IE, JE, KE)
    hx = zeros(Float64, IE, JE, KE)
    hy = zeros(Float64, IE, JE, KE)
    hz = zeros(Float64, IE, JE, KE)
    gax = ones(Float64, IE, JE, KE)
    gay = ones(Float64, IE, JE, KE)
    gaz = ones(Float64, IE, JE, KE)

    # Specify the dipole
    for k in (kc - 10):(kc + 10)
        if k >= 1 && k <= KE
            gaz[ic, jc, k] = 0.0
        end
    end
    if kc >= 1 && kc <= KE
        gaz[ic, jc, kc] = 1.0
    end

    # Prepare the plot
    should_display = true  # Set to false if running in a non-interactive environment

    # Create grid coordinates for plotting
    x = (1:IE) .* ddx
    y = (1:JE) .* ddx

    # Create coordinate matrices
    X = repeat(x', JE, 1)
    Y = repeat(y, 1, IE)

    # Time-stepping loop with animation
    anim = @animate for n = 1:nsteps
        T = n  # Time variable matching the C code

        # --- Update Dx field ---
        for i = 2:IE
            for j = 2:JE
                for k = 2:KE
                    dx[i, j, k] += 0.5 * (hz[i, j, k] - hz[i, j - 1, k] - hy[i, j, k] + hy[i, j, k - 1])
                end
            end
        end

        # --- Update Dy field ---
        for i = 2:IE
            for j = 2:JE
                for k = 2:KE
                    dy[i, j, k] += 0.5 * (hx[i, j, k] - hx[i, j, k - 1] - hz[i, j, k] + hz[i - 1, j, k])
                end
            end
        end

        # --- Update Dz field ---
        for i = 2:IE
            for j = 2:JE
                for k = 2:KE
                    dz[i, j, k] += 0.5 * (hy[i, j, k] - hy[i - 1, j, k] - hx[i, j, k] + hx[i, j - 1, k])
                end
            end
        end

        # --- Add the source at the gap ---
        pulse = exp(-0.5 * ((t0 - T) / spread)^2)
        if ic >= 1 && ic <= IE && jc >= 1 && jc <= JE && kc >= 1 && kc <= KE
            dz[ic, jc, kc] = pulse  # Assign the pulse directly
        end

        # --- Calculate the E from D field ---
        for i = 2:(IE - 1)
            for j = 2:(JE - 1)
                for k = 2:(KE - 1)
                    ex[i, j, k] = gax[i, j, k] * dx[i, j, k]
                    ey[i, j, k] = gay[i, j, k] * dy[i, j, k]
                    ez[i, j, k] = gaz[i, j, k] * dz[i, j, k]
                end
            end
        end

        # Print time step and Ez at dipole
        println("Time step: $T, Ez at dipole: $(ez[ic, jc, kc])")

        # --- Update Hx field ---
        for i = 2:IE
            for j = 2:(JE - 1)
                for k = 2:(KE - 1)
                    hx[i, j, k] += 0.5 * (ey[i, j, k + 1] - ey[i, j, k] - ez[i, j + 1, k] + ez[i, j, k])
                end
            end
        end

        # --- Update Hy field ---
        for i = 2:(IE - 1)
            for j = 2:JE
                for k = 2:(KE - 1)
                    hy[i, j, k] += 0.5 * (ez[i + 1, j, k] - ez[i, j, k] - ex[i, j, k + 1] + ex[i, j, k])
                end
            end
        end

        # --- Update Hz field ---
        for i = 2:(IE - 1)
            for j = 2:(JE - 1)
                for k = 2:KE
                    hz[i, j, k] += 0.5 * (ex[i, j + 1, k] - ex[i, j, k] - ey[i + 1, j, k] + ey[i, j, k])
                end
            end
        end

        # --- Visualization ---
        if n % 2 == 0 || n == 1
            # Extract Ez at z = kc
            Ez_slice = ez[:, :, kc]

            # Create a 3D surface plot
            surface(
                X, Y, Ez_slice',
                clims = (-0.005, 0.005),  # Adjust clims based on field values
                title = "Ez Field at z = $kc, Time Step $n",
                xlabel = "X (m)",
                ylabel = "Y (m)",
                zlabel = "Ez (V/m)",
                c = :RdBu,
                colorbar = true,
                camera = (30, 30)  # Adjust viewing angle if desired
            )
        end
    end  # End of @animate

    # Display or save the animation
    if save_gif
        # Save the animation as a GIF file
        gif(anim, "fdtd_3d_dipole.gif", fps = 10)
        println("Animation saved as fdtd_3d_dipole.gif")
    else
        # Display the animation
        if should_display
            display(anim)
        end
    end

    println("Simulation complete. Final time step: $nsteps")
end


function fdtd_3d_sphere_pml(nsteps::Int; save_gif::Bool = false)
    """
    3D FDTD Simulation with a dielectric sphere in Julia with PML

    """
    # Constants
    IE = 100
    JE = 100
    KE = 40
    ia = 7
    ja = 7
    ka = 7
    NFREQS = 3

    # Physical constants
    epsz = 8.8e-12  # Vacuum permittivity (F/m)
    muz = 4π * 1e-7  # Vacuum permeability (H/m)
    pi = π
    ddx = 0.01  # Cell size (m)
    dt = ddx / (6e8)  # Time step size (s)

    # Center indices
    ic = div(IE, 2)
    jc = div(JE, 2)
    kc = div(KE, 2)
    ib = IE - ia - 1
    jb = JE - ja - 1
    kb = KE - ka - 1

    # Arrays for fields
    dx = zeros(Float64, IE, JE, KE)
    dy = zeros(Float64, IE, JE, KE)
    dz = zeros(Float64, IE, JE, KE)
    ex = zeros(Float64, IE, JE, KE)
    ey = zeros(Float64, IE, JE, KE)
    ez = zeros(Float64, IE, JE, KE)
    hx = zeros(Float64, IE, JE, KE)
    hy = zeros(Float64, IE, JE, KE)
    hz = zeros(Float64, IE, JE, KE)
    ix = zeros(Float64, IE, JE, KE)
    iy = zeros(Float64, IE, JE, KE)
    iz = zeros(Float64, IE, JE, KE)

    # Arrays for PML, Fourier transforms, and more
    gax = ones(Float64, IE, JE, KE)
    gay = ones(Float64, IE, JE, KE)
    gaz = ones(Float64, IE, JE, KE)
    gbx = zeros(Float64, IE, JE, KE)
    gby = zeros(Float64, IE, JE, KE)
    gbz = zeros(Float64, IE, JE, KE)
    idx = zeros(Float64, IE, JE, KE)
    idy = zeros(Float64, IE, JE, KE)
    idz = zeros(Float64, IE, JE, KE)
    ihx = zeros(Float64, IE, JE, KE)
    ihy = zeros(Float64, IE, JE, KE)
    ihz = zeros(Float64, IE, JE, KE)

    # Incident fields
    ez_inc = zeros(Float64, JE)
    hx_inc = zeros(Float64, JE)
    ez_low_m1 = 0.0
    ez_low_m2 = 0.0
    ez_high_m1 = 0.0
    ez_high_m2 = 0.0

    # Fourier Transform arrays
    real_in = zeros(Float64, NFREQS)
    imag_in = zeros(Float64, NFREQS)
    real_pt = zeros(Float64, NFREQS, IE, JE)
    imag_pt = zeros(Float64, NFREQS, IE, JE)
    amp = zeros(Float64, IE, JE)
    phase = zeros(Float64, IE, JE)
    freq = [10e6, 100e6, 433e6]
    arg = 2 * pi * freq * dt

    # Boundary condition arrays
    gi1 = zeros(Float64, IE)
    gi2 = ones(Float64, IE)
    gi3 = ones(Float64, IE)
    gj1 = zeros(Float64, JE)
    gj2 = ones(Float64, JE)
    gj3 = ones(Float64, JE)
    gk1 = zeros(Float64, KE)
    gk2 = ones(Float64, KE)
    gk3 = ones(Float64, KE)
    fi1 = zeros(Float64, IE)
    fi2 = ones(Float64, IE)
    fi3 = ones(Float64, IE)
    fj1 = zeros(Float64, JE)
    fj2 = ones(Float64, JE)
    fj3 = ones(Float64, JE)
    fk1 = zeros(Float64, KE)
    fk2 = ones(Float64, KE)
    fk3 = ones(Float64, KE)

    # Sphere parameters
    numsph = 1
    radius = [0.0, 5.0]  # Indexing from 1, radius[1] is the sphere radius
    epsilon = [1.0, 3.0]  # epsilon[1] is the sphere permittivity
    sigma = [0.0, 0.0]    # sigma[1] is the sphere conductivity

    # Set up gax, gbx, gay, gby, gaz, gbz based on the sphere
    for i in ia:ib
        for j in ja:jb
            for k in ka:kb
                # For ex
                xdist = (ic - i - 0.5)
                ydist = (jc - j)
                zdist = (kc - k)
                dist = sqrt(xdist^2 + ydist^2 + zdist^2)
                eps = epsilon[1]
                cond = sigma[1]
                if dist <= radius[2]
                    eps = epsilon[2]
                    cond = sigma[2]
                end
                gax[i, j, k] = 1.0 / (eps + (cond * dt / epsz))
                gbx[i, j, k] = cond * dt / epsz
                # For ey
                xdist = (ic - i)
                ydist = (jc - j - 0.5)
                zdist = (kc - k)
                dist = sqrt(xdist^2 + ydist^2 + zdist^2)
                eps = epsilon[1]
                cond = sigma[1]
                if dist <= radius[2]
                    eps = epsilon[2]
                    cond = sigma[2]
                end
                gay[i, j, k] = 1.0 / (eps + (cond * dt / epsz))
                gby[i, j, k] = cond * dt / epsz
                # For ez
                xdist = (ic - i)
                ydist = (jc - j)
                zdist = (kc - k - 0.5)
                dist = sqrt(xdist^2 + ydist^2 + zdist^2)
                eps = epsilon[1]
                cond = sigma[1]
                if dist <= radius[2]
                    eps = epsilon[2]
                    cond = sigma[2]
                end
                gaz[i, j, k] = 1.0 / (eps + (cond * dt / epsz))
                gbz[i, j, k] = cond * dt / epsz
            end
        end
    end

    # Time-stepping loop
    t0 = 40.0
    spread = 10.0
    T = 0.0  # Initialize T before the loop

    # Initialize animation
    anim = Animation()

    # For animation
    for n in 1:nsteps
        T += 1.0  # Increment T

        # --- Incident buffer calculation ---
        for j in 2:JE-1
            ez_inc[j] += 0.5 * (hx_inc[j-1] - hx_inc[j])
        end

        # --- Source pulse ---
        pulse = exp(-0.5 * ((t0 - T) / spread)^2)
        ez_inc[4] = pulse  # Adjusted index for Julia (ez_inc[3] in C, but Julia is 1-based)

        # --- Boundary conditions for incident buffer ---
        ez_inc[1] = ez_low_m2
        ez_low_m2 = ez_low_m1
        ez_low_m1 = ez_inc[2]

        ez_inc[JE] = ez_high_m2
        ez_high_m2 = ez_high_m1
        ez_high_m1 = ez_inc[JE-1]

        # --- Update Dx field ---
        for i in 2:IE-1
            for j in 2:JE-1
                for k in 2:KE-1
                    curl_h = hz[i, j, k] - hz[i, j-1, k] - hy[i, j, k] + hy[i, j, k-1]
                    idx[i, j, k] += curl_h
                    dx[i, j, k] = gj3[j] * gk3[k] * dx[i, j, k] + gj2[j] * gk2[k] * 0.5 * (curl_h + gi1[i] * idx[i, j, k])
                end
            end
        end

        # --- Update Dy field ---
        for i in 2:IE-1
            for j in 2:JE-1
                for k in 2:KE-1
                    curl_h = hx[i, j, k] - hx[i, j, k-1] - hz[i, j, k] + hz[i-1, j, k]
                    idy[i, j, k] += curl_h
                    dy[i, j, k] = gi3[i] * gk3[k] * dy[i, j, k] + gi2[i] * gk2[k] * 0.5 * (curl_h + gj1[j] * idy[i, j, k])
                end
            end
        end

        # --- Incident Dy ---
        for i in ia:ib
            for j in ja:jb-1
                dy[i, j, ka] -= 0.5 * hx_inc[j]
                dy[i, j, kb+1] += 0.5 * hx_inc[j]
            end
        end

        # --- Update Dz field ---
        for i in 2:IE-1
            for j in 2:JE-1
                for k in 1:KE-1
                    curl_h = hy[i, j, k] - hy[i-1, j, k] - hx[i, j, k] + hx[i, j-1, k]
                    idz[i, j, k] += curl_h
                    dz[i, j, k] = gi3[i] * gj3[j] * dz[i, j, k] + gi2[i] * gj2[j] * 0.5 * (curl_h + gk1[k] * idz[i, j, k])
                end
            end
        end

        # --- Incident Dz ---
        for i in ia:ib
            for k in ka:kb
                dz[i, ja, k] += 0.5 * hx_inc[ja]
                dz[i, jb, k] -= 0.5 * hx_inc[jb]
            end
        end

        # --- Update E fields from D fields ---
        for i in 2:IE-1
            for j in 2:JE-1
                for k in 2:KE-1
                    ex[i, j, k] = gax[i, j, k] * (dx[i, j, k] - ix[i, j, k])
                    ix[i, j, k] += gbx[i, j, k] * ex[i, j, k]
                    ey[i, j, k] = gay[i, j, k] * (dy[i, j, k] - iy[i, j, k])
                    iy[i, j, k] += gby[i, j, k] * ey[i, j, k]
                    ez[i, j, k] = gaz[i, j, k] * (dz[i, j, k] - iz[i, j, k])
                    iz[i, j, k] += gbz[i, j, k] * ez[i, j, k]
                end
            end
        end

        # --- Fourier Transform of Ez ---
        for m in 1:NFREQS
            for j in 1:JE
                for i in 1:IE
                    real_pt[m, i, j] += cos(arg[m] * T) * ez[i, j, kc]
                    imag_pt[m, i, j] += sin(arg[m] * T) * ez[i, j, kc]
                end
            end
        end

        # --- Update incident Hx field ---
        for j in 1:JE-1
            hx_inc[j] += 0.5 * (ez_inc[j] - ez_inc[j+1])
        end

        # --- Update Hx field ---
        for i in 1:IE-1
            for j in 1:JE-1
                for k in 1:KE-1
                    curl_e = ey[i, j, k+1] - ey[i, j, k] - ez[i, j+1, k] + ez[i, j, k]
                    ihx[i, j, k] += curl_e
                    hx[i, j, k] = fj3[j] * fk3[k] * hx[i, j, k] + fj2[j] * fk2[k] * 0.5 * (curl_e + fi1[i] * ihx[i, j, k])
                end
            end
        end

        # --- Incident Hx ---
        for i in ia:ib
            for k in ka:kb
                hx[i, ja-1, k] += 0.5 * ez_inc[ja]
                hx[i, jb, k] -= 0.5 * ez_inc[jb]
            end
        end

        # --- Update Hy field ---
        for i in 1:IE-1
            for j in 1:JE-1
                for k in 1:KE-1
                    curl_e = ez[i+1, j, k] - ez[i, j, k] - ex[i, j, k+1] + ex[i, j, k]
                    ihy[i, j, k] += curl_e
                    hy[i, j, k] = fi3[i] * fk3[k] * hy[i, j, k] + fi2[i] * fk2[k] * 0.5 * (curl_e + fj1[j] * ihy[i, j, k])
                end
            end
        end

        # --- Incident Hy ---
        for j in ja:jb
            for k in ka:kb
                hy[ia-1, j, k] -= 0.5 * ez_inc[j]
                hy[ib, j, k] += 0.5 * ez_inc[j]
            end
        end

        # --- Update Hz field ---
        for i in 1:IE-1
            for j in 1:JE-1
                for k in 1:KE
                    curl_e = ex[i, j+1, k] - ex[i, j, k] - ey[i+1, j, k] + ey[i, j, k]
                    ihz[i, j, k] += curl_e
                    hz[i, j, k] = fi3[i] * fj3[j] * hz[i, j, k] + fi2[i] * fj2[j] * 0.5 * (curl_e + fk1[k] * ihz[i, j, k])
                end
            end
        end

        # --- 3D Surface Visualization ---
        if n % 2 == 0 || n == 1
            # Create X, Y, Z grids
            x = 1:IE
            y = 1:JE
            z = ez[:, :, kc]  # Get the Ez field at z = kc (middle plane)

            # Plot the 3D surface
            p = surface(x, y, z', xlabel = "X", ylabel = "Y", zlabel = "Ez",
                        title = "Ez Field at z = $kc, Time Step $n", c = :RdBu,zlims=(-0.3, 1.2))
            frame(anim, p)
        end
    end  # End of time-stepping loop

    # Save or display the animation
    if save_gif
        gif(anim, "fdtd_3d_sphere.gif", fps = 10)
        println("Animation saved as fdtd_3d_sphere.gif")
    else
        display(anim)
    end

    println("Simulation complete. Final time step: $nsteps")
end



end 