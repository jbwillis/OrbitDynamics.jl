# Functions and structs for simulating orbits

struct DynamicsParameters
    # fixed parameters
	mu::Float64
    J2::Float64
    R_earth::Float64
    m_satellite::Float64
    rho::Float64
    C_D::Float64
    A::Float64
	distance_scale::Float64
	time_scale::Float64
    
    function DynamicsParameters(;
            m_satellite_si=1.0, # kg
            rho_si=1e-12, # kg/m^3
            C_D=1.0, # nondimensional
            A_si=(.3*.1), # m^2
			mu_earth_si=3.986004415e14, # m^3/s^2
			J2=0.0010826358191967, # nondimensional
			R_earth_si=6.378136300e6, # m
			distance_scale=1.0, # m/distance_unit
			time_scale=1.0, # s/time_unit
			)
        return new(
		    mu_earth_si * (time_scale^2) / (distance_scale^3), # mu
            J2,
            R_earth_si / distance_scale,
            m_satellite_si,
			rho_si * (distance_scale^3),
            C_D,
			A_si / (distance_scale^2),
			distance_scale,
			time_scale
			)
    end
end

"""
Taken from Spacecraft Attitude Determination and Control - works well
"""
function gravity_force_ECI(x, p::DynamicsParameters)
	r = x[1:3]
    r_mag = sqrt(r'*r)
    a_g = -(p.mu/r_mag^3) * r
	a_J2 = gravity_perturbation_ECI(x, p)
    F_g = p.m_satellite * (a_g + a_J2)

    return F_g
end

function gravity_perturbation_ECI(x, p::DynamicsParameters)
	r = x[1:3]
    r_mag = sqrt(r'*r)
    a_J2 = -(3.0/2.0) * p.J2 * (p.mu/r_mag^2) * (p.R_earth/r_mag)^2 * [ 
        (1 - 5*(r[3]/r_mag)^2) * r[1]/r_mag,        
        (1 - 5*(r[3]/r_mag)^2) * r[2]/r_mag,        
        (3 - 5*(r[3]/r_mag)^2) * r[3]/r_mag]

	return a_J2
end
    

function orbit_dynamics_ECI_state(x, p::DynamicsParameters, t)
    
    r = x[1:3]
    r_mag = sqrt(r'*r)
    v = x[4:6]
    v_mag = sqrt(v'*v)            
        
    # Gravity        
    F_g = gravity_force_ECI(x, p)    
    
    # Drag
    F_d = -0.5 * p.rho * p.C_D * p.A * v_mag .* v
    
	x_dot = [v; (1.0/p.m_satellite) .* (F_g .+ F_d)]
   
	return x_dot
end

function orbit_dynamics_ECI_drag_control(x, p::DynamicsParameters, t)
    r = x[1:3]
    v = x[4:6]
    u = x[7]
    
    r_mag = sqrt(r'*r)
    v_mag = sqrt(v'*v)            
        
    # Gravity        
    F_g = gravity_force_ECI(x, p)    
    
    # Drag
    F_d = -0.5 * p.rho * p.C_D * (u * p.A) * v_mag .* v
    
    p_dot = v
    v_dot = (1.0/p.m_satellite) .* (F_g .+ F_d)
    u_dot = 0.0
	x_dot = [p_dot; v_dot; u_dot]
   
	return x_dot
end


function solve_orbit_dynamics_ECI_state(x0, dp::DynamicsParameters, t_end, h=0.1)
    
	x, t = integrate_RK4(orbit_dynamics_ECI_state, x0, dp, 0.0, t_end, h)

    return x, t
end

function circular_orbit_initial_conditions(orbit_altitude, inclination, p::DynamicsParameters)
    
    r_orbit = orbit_altitude + p.R_earth
    
    v_orbit = sqrt(p.mu/r_orbit)
    
    r_vec = [r_orbit, 0.0, 0.0] # start at equator
    v_vec = [0.0, v_orbit * cos(inclination), v_orbit * sin(inclination)] # velocity along path of orbit
    
    x = vcat(r_vec, v_vec)
    
    return x
end

"""
Use the modified equinoctial elements to simulate.
Allows large time steps and prevents singularities for e=0, i=0,90 degrees

See: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf

x = [p, f, g, h, k, L]
"""
function orbit_dynamics_equinoctial(x, dp::DynamicsParameters, t)

	p, f, g, h, k, L = x

	u_J2 = gravity_perturbation_equinoctial(x, dp)
	u_drag = drag_perturbation_equinoctial(x, dp)

	u_R, u_T, u_N = (u_J2 + u_drag)

	w = 1 + (f * cos(L)) + (g * sin(L))
	s_sq = 1 + h^2 + k^2

	sqrtpu = sqrt(p/dp.mu)

	p_dot = 2 * (p / w) * sqrtpu * u_T

	f_dot = sqrtpu * (( sin(L) * u_R) + ((1 / w) * ((w + 1) * cos(L) + f) * u_T) - ((g / w) * (h * sin(L) - k * cos(L)) * u_N))

	g_dot = sqrtpu * ((-cos(L) * u_R) + ((1 / w) * ((w + 1) * sin(L) + g) * u_T) + ((f / w) * (h * sin(L) - k * cos(L)) * u_N))

	h_dot = sqrtpu * (s_sq / (2 * w)) * cos(L) * u_N

	k_dot = sqrtpu * (s_sq / (2 * w)) * sin(L) * u_N

	L_dot = ( sqrt(dp.mu * p) * (w / p)^2) + ((1 / w) * sqrtpu * (h * sin(L) - k * cos(L)) * u_N)

	x_dot = [p_dot, f_dot, g_dot, h_dot, k_dot, L_dot]

	return x_dot

end

function solve_orbit_dynamics_equinoctial(x0, dp::DynamicsParameters, t_end, h=1.0)
    
	x, t = integrate_RK4(orbit_dynamics_equinoctial, x0, dp, 0.0, t_end, h)

    return x, t
end

"""
J2 perturbing gravity forces computed using equinoctial coordinates
Source: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
"""
function gravity_perturbation_equinoctial(x, dp)

	p, f, g, h, k, L = x

	w = 1 + f * cos(L) + g * sin(L)
	r = p / w

	u_J2_R = ((-3  * dp.mu * dp.J2 * dp.R_earth^2) / (2 * r^4)) * (1 - 12 * (h * sin(L) - k * cos(L))^2 / (1 + h^2 + k^2)^2)
	u_J2_T = ((-12 * dp.mu * dp.J2 * dp.R_earth^2) / (r^4))     * ((h * sin(L) - k * cos(L))*(h * cos(L) + k * sin(L)) / (1 + h^2 + k^2)^2)
	u_J2_N = ((-6  * dp.mu * dp.J2 * dp.R_earth^2) / (r^4))     * ((1 - h^2 - k^2)*(h * sin(L) - k * cos(L)) / (1 + h^2 + k^2)^2)

	u_J2 = [u_J2_R, u_J2_T, u_J2_N]
	return u_J2
end

"""
Aerodynamic drag accelerations using equinoctial coordinates
"""
function drag_perturbation_equinoctial(x, dp::DynamicsParameters)

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

function orbit_dynamics_classical_elements(x, dp::DynamicsParameters, t)
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
	omega_dot = (1/(e*h)) * ((-p * cos(theta) * u_R) + ((p + r) * sin(theta) * u_T)) - ((r * sin(theta + omega) * cos(i)) * u_N / (h * sin(i)))
	Omega_dot = (r * sin(theta + omega) * u_N) / (h * sin(i))
	theta_dot = (h/r^2) + (1/(e*h)) * ((p * cos(theta) * u_R) - ((p + r) * sin(theta) * u_T))

	x_dot = [a_dot, e_dot, i_dot, omega_dot, Omega_dot, theta_dot]

	return x_dot
end

function solve_orbit_dynamics_classical_elements(x0, dp::DynamicsParameters, t_end, h=1.0)
    
	x, t = integrate_RK4(orbit_dynamics_classical_elements, x0, dp, 0.0, t_end, h)

    return x, t
end

"""
Source: Orbital Mechanics for Engineering Students, pg 178
"""
function gravity_perturbation_classical(x, dp)
	a, e, i, omega, Omega, theta = x
	p = a*(1-e^2)
	r = p / (1 + e*cos(theta))
	u = omega + theta

	c_J2 = (-3 * dp.mu * dp.J2 * dp.R_earth^2)/(2 * r^4)

	u_J2_R = c_J2 * (1 - 3 * sin(i)^2 * sin(u)^2)
	u_J2_T = 2*c_J2 * (sin(i)^2 * sin(u) * cos(u))
	u_J2_N = 2*c_J2 * (sin(i) * cos(i) * sin(u))

	u_J2 = [u_J2_R, u_J2_T, u_J2_N]

	return u_J2
end

"""
Source: Position and velocity perturbations in the orbital frame in terms of classical element perturbations
"""
function drag_perturbation_classical(x, dp)

	a, e, i, omega, Omega, theta = x

	n = sqrt(dp.mu/(a^3)) # mean motion

	v_R = (n * a * e * sin(theta))/sqrt(1 - e^2)
	v_T = (n * a * (1 + e*cos(theta)))/sqrt(1 - e^2)

	return drag_perturbation_RTN(v_R, v_T, dp)
end

"""
Use Runge-Kutta 4 to integrate the dynamics from `t_start` to `t_end` 
using `x0` as the initial condition, 
`p` as parameters for the dynamics function,
and `h` as the fixed timestep.
"""
function integrate_RK4(dynamics, x0, p, t_start, t_end, h)

	N = Int(round((t_end - t_start) / h) + 1)

	# allocate arrays
	x = zeros(length(x0), N)
	t = zeros(N)

	# initial condition
	x[:,1] .= x0
	t[1] = 0.0

	# integrate
	for i=2:N
		t[i] = h*(i-1)
		x[:,i] = step_RK4(dynamics, x[:,i-1], p, t[i], h)
	end

	return x, t
end

"""
Use Runge-Kutta 4 to integrate the dynamics at `(x, t)`, using parameters `p` and timestep `h`.
"""
function step_RK4(dynamics, x_k, p, t, h)

	k1 = dynamics(x_k, p, t)
	k2 = dynamics(x_k .+ 0.5 .* h .* k1, p, t + 0.5 * h)
	k3 = dynamics(x_k .+ 0.5 .* h .* k2, p, t + 0.5 * h)
	k4 = dynamics(x_k .+ h .* k3, p, t + h)

	x_kp1 = x_k .+ (1.0/6.0) * h .* (k1 .+ (2 .* k2) .+ (2 .* k3) .+ k4)

	return x_kp1

end

