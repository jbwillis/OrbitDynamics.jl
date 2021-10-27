"""
Functions and structs for simulating orbits
"""

struct DynamicsParameters
    # fixed parameters
    G_gravity
    M_earth
	mu
    J2
    R_earth
    m_satellite
    rho
    C_D
    A
    
    function DynamicsParameters(;
            m_satellite=1.0,
            rho=1e-12,
            C_D=1.0, 
            A=(.3*.1), # side of a 3U
            G_gravity=6.674e-11,
            M_earth=5.9722e24, # kg
            J2=0.1082626925638815e-2, 
            R_earth=6378137.)
        return new(
            G_gravity,
            M_earth,
			mu = G_gravity * M_earth,
            J2,
            R_earth,
            m_satellite,
            rho,
            C_D,
            A)
    end
end

function gravity_force_from_r(r, p::DynamicsParameters)
    """
    Taken from Spacecraft Attitude Determination and Control - works well
    """
    r_mag = sqrt(r'*r)
    a_g = -(p.mu/r_mag^3) * r
    a_J2 = -(3.0/2.0) * p.J2 * (p.mu/r_mag^2) * (p.R_earth/r_mag)^2 * [ 
        (1 - 5*(r[3]/r_mag)^2) * r[1]/r_mag,        
        (1 - 5*(r[3]/r_mag)^2) * r[2]/r_mag,        
        (3 - 5*(r[3]/r_mag)^2) * r[3]/r_mag]
    
    F_g = p.m_satellite * (a_g + a_J2)

    return F_g
end
    

function orbit_dynamics_ECI_state!(x_dot, x, p::DynamicsParameters, t)
    
    r = x[1:3]
    r_mag = sqrt(r'*r)
    v = x[4:6]
    v_mag = sqrt(v'*v)            
        
    # Gravity        
    F_g = gravity_force_from_r(r, p)    
    
    # Drag
    F_d = -0.5 * p.rho * p.C_D * p.A * v_mag .* v
    
    x_dot[1:3] = v
    x_dot[4:6] = (1.0/p.m_satellite) .* (F_g .+ F_d)
    
end

function solve_orbit_dynamics_ECI_state(x0, p::DynamicsParameters, t_end)
    
    t_span = (0.0, t_end)    
    
    prob = ODEProblem(orbit_dynamics!, x0, t_span, p)
    sol = solve(prob, saveat=t_end/1e5, abstol=1e-8, reltol=1e-8);
    
    r_sol = sol[1:3,:]
    v_sol = sol[4:6,:]
    t_sol = sol.t
    return r_sol, v_sol, t_sol
end

function circular_orbit_initial_conditions(orbit_altitude, inclination, p::DynamicsParameters)
    
    r_orbit = orbit_altitude + p.R_earth
    
    v_orbit = sqrt(p.mu/r_orbit)
    
    r_vec = [r_orbit, 0.0, 0.0] # start at equator
    v_vec = [0.0, v_orbit * cos(inclination), v_orbit * sin(inclination)] # velocity along path of orbit
    
    x = vcat(r_vec, v_vec)
    
    return x
end

function orbit_dynamics_equinoctial!(x_dot, x, dp::DynamicsParameters, t)
"""
Use the modified equinoctial elements to simulate.
Allows large time steps and prevents singularities for e=0, i=0,90 degrees

See: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf

x = [p, f, g, h, k, L]
"""

p, f, g, h, k, L = x

u_J2 = gravity_perturbation_equinoctial(x, dp)
u_drag = drag_perturbation_equinoctial(x, dp)

u_R, u_T, u_N = (u_J2 + u_drag)

w = 1 + f * cos(L) + g * sin(L)
s_sq = 1 + h^2 + k^2

p_dot = 2 * (p / w) * sqrt(p / dp.mu) * u_R
f_dot = ((sqrt(p / dp.mu) * sin(L) * u_R)
		 + (sqrt(p / dp.mu) * (1 / w) * ((w + 1) * cos(L) + f) * u_T)
		 - (sqrt(p / dp.mu) * (g / w) * (h * sin(L) - k * cos(L)) * u_N))
g_dot = ((-sqrt(p / dp.mu) * cos(L) * u_R)
		 + (sqrt(p / dp.mu) * ((w + 1) * sin(L) + g) * u_T)
		 - (sqrt(p / dp.mu) * (f / w) * (h * sin(L) - k * cos(L)) * u_N))
h_dot = sqrt(p/dp.mu) * (s_sq / (2 * w)) * cos(L) * u_N
k_dot = sqrt(p/dp.mu) * (s_sq / (2 * w)) * sin(L) * u_N
L_dot = ((sqrt(dp.mu * p) * (w / p)^2)
		 + (sqrt(p / dp.mu) * (h * sin(L) - k * cos(L)) * u_N))

x_dot[1] = p_dot
x_dot[2] = f_dot
x_dot[3] = g_dot
x_dot[4] = h_dot
x_dot[5] = k_dot
x_dot[6] = L_dot

end

function gravity_perturbation_equinoctial(x, dp)
	"""
	J2 perturbing gravity forces computed using equinoctial coordinates
	"""

	p, f, g, h, k, L = x

	w = 1 + f * cos(L) + g * sin(L)
	r = p / w

	u_J2_R = ((-3 * dp.mu * dp.J2 * dp.R_earth^2) / (2 * r^4)) * (1 - 12 * (h * sin(L) - k * cos(L))^2 / (1 + h^2 + k^2)^2)
	u_J2_T = ((-12 * dp.mu * dp.J2 * dp.R_earth^2) / (2 * r^4)) * ((h * sin(L) - k * cos(L))*(h * cos(L) + k * sin(L)) / (1 + h^2 + k^2)^2)
	u_J2_N = ((-6 * dp.mu * dp.J2 * dp.R_earth^2) / (2 * r^4)) * ((1 - h^2 - k^2)*(h * sin(L) - k * cos(L)) / (1 + h^2 + k^2)^2)

	u_J2 = [u_J2_R, u_J2_T, u_J2_N]
	return u_J2
end

function drag_perturbation_equinoctial(x, dp)
	"""
	Aerodynamic drag accelerations using equinoctial coordinates
	"""

	p, f, g, h, k, L = x

	v_R = sqrt(dp.mu / p) * (f*sin(L) - g*cos(L))
	v_T = sqrt(dp.mu / p) * (1 + f*cos(L) + g*sin(L))
	
	v = sqrt(v_R^2 + v_T^2)

	u_drag_R = -0.5 * dp.rho * dp.A * dp.C_D * v * v_R
	u_drag_T = -0.5 * dp.rho * dp.A * dp.C_D * v * v_T
	u_drag_N = 0.0

	u_drag = [u_drag_R, u_drag_T, u_drag_N]

	return u_drag
end
