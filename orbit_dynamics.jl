# Functions and structs for simulating orbits

using DifferentialEquations

struct DynamicsParameters
    # fixed parameters
    G_gravity::Float64
    M_earth::Float64
	mu::Float64
    J2::Float64
    R_earth::Float64
    m_satellite::Float64
    rho::Float64
    C_D::Float64
    A::Float64
    
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
			G_gravity * M_earth,
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
    
    prob = ODEProblem(orbit_dynamics_ECI_state!, x0, t_span, p)
    sol = solve(prob, saveat=t_end/1e4, abstol=1e-8, reltol=1e-8);
    
    x_sol = sol[1:6,:]
    t_sol = sol.t
    return x_sol, t_sol
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

function solve_orbit_dynamics_equinoctial(x0, dp::DynamicsParameters, t_end)
    
    t_span = (0.0, t_end)    
    
    prob = ODEProblem(orbit_dynamics_equinoctial!, x0, t_span, dp)
    sol = solve(prob, saveat=t_end/1e4, abstol=1e-8, reltol=1e-8);
    
    x_sol = sol[1:6,:]
    t_sol = sol.t
    return x_sol, t_sol
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

function drag_perturbation_equinoctial(x, dp::DynamicsParameters)
	"""
	Aerodynamic drag accelerations using equinoctial coordinates
	"""

	p, f, g, h, k, L = x

	v_R = sqrt(dp.mu / p) * (f*sin(L) - g*cos(L))
	v_T = sqrt(dp.mu / p) * (1 + f*cos(L) + g*sin(L))
	
	return drag_perturbation_RTN(v_R, v_T, dp)
end

function drag_perturbation_RTN(v_R, v_T, dp::DynamicsParameters)
	v = sqrt(v_R^2 + v_T^2)

	u_drag_R = -0.5 * dp.rho * dp.A * dp.C_D * v * v_R
	u_drag_T = -0.5 * dp.rho * dp.A * dp.C_D * v * v_T
	u_drag_N = 0.0

	u_drag = [u_drag_R, u_drag_T, u_drag_N]

	return u_drag
end

function orbit_dynamics_classical_elements!(x_dot, x, dp::DynamicsParameters, t)
	a, e, i, omega, Omega, theta = x

	p = a*(1-e^2)
	r = p / (1 + e*cos(theta))
	h = sqrt(dp.mu * p)

	u_J2 = gravity_perturbation_classical(x, dp)
	u_drag = drag_perturbation_classical(x, dp)

	u_R, u_T, u_N = (u_J2 + u_drag)

	a_dot = 2 * (a^2 / h) * ((e*sin(theta)*u_R) + ((p/r) * u_T))
	e_dot = (1/h) * (p*sin(theta)*u_R + (((p+r)*cos(theta) + r*e)*u_T))
	i_dot = r * cos(theta + omega) * u_N / h
	omega_dot = (1/e*h) * ((-p * cos(theta) * u_R) + ((p + r) * sin(theta) * u_T)) - ((r * sin(theta + omega) * cos(i)) * u_N / (h * sin(i)))
	Omega_dot = (r * sin(theta + omega) * u_N) / (h * sin(i))
	theta_dot = (h/r^2) + (1/(e*h)) * ((p * cos(theta) * u_R) - ((p + r) * sin(theta) * u_T))

	x_dot[1] = a_dot
	x_dot[2] = e_dot
	x_dot[3] = i_dot
	x_dot[4] = omega_dot
	x_dot[5] = Omega_dot
	x_dot[6] = theta_dot

end

function solve_orbit_dynamics_classical_elements(x0, dp::DynamicsParameters, t_end)
    
    t_span = (0.0, t_end)    
    
    prob = ODEProblem(orbit_dynamics_classical_elements!, x0, t_span, dp)
    sol = solve(prob, saveat=t_end/1e4, abstol=1e-8, reltol=1e-8);
    
    x_sol = sol[1:6,:]
    t_sol = sol.t
    return x_sol, t_sol
end

function gravity_perturbation_classical(x, dp)
	a, e, i, omega, Omega, theta = x
	p = a*(1-e^2)
	r = p / (1 + e*cos(theta))
	u = omega + theta

	c_J2 = (-3 * dp.mu * dp.J2 * dp.R_earth^2)/(2 * r^4)

	u_J2_R = c_J2 * (1 - 3 * sin(i)^2 * sin(u)^2)
	u_J2_T = c_J2 * (sin(i)^2 * sin(u) * cos(u))
	u_J2_N = c_J2 * (sin(i) * cos(i) * sin(u))

	u_J2 = [u_J2_R, u_J2_T, u_J2_N]

	return u_J2
end

function drag_perturbation_classical(x, dp)
	"""
	Source: Position and velocity perturbations in the orbital frame in terms of classical element perturbations 
	"""

	a, e, i, omega, Omega, theta = x

	n = sqrt(dp.mu/(a^3)) # mean motion

	v_R = (n * a * e * sin(theta))/sqrt(1 - e^2)
	v_T = (n * a * (1 + e*cos(theta)))/sqrt(1 - e^2)

	return drag_perturbation_RTN(v_R, v_T, dp)
end
