# Code for converting between different orbit representations

using LinearAlgebra

include("orbit_dynamics.jl")

function equinoctial_to_classical_elements(x_eq)
	"""
	Reference: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
	"""
	p, f, g, h, k, L = x_eq

	a = p / (1 - f^2 - g^2)
	e = sqrt(f^2 + g^2)
	i = atan(2*sqrt(h^2 + k^2), 1 - h^2 - k^2)
	omega = atan(g*h - f*k, f*h + g*k)
	Omega = atan(k, h)
	theta = L - (Omega + omega)

	omega = mod2pi(omega)
	Omega = mod2pi(Omega)
	theta = mod2pi(theta)

	return [a, e, i, omega, Omega, theta]
end

function classical_to_equinoctial_elements(x_cl)
	"""
	Reference: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
	"""
	a, e, i, omega, Omega, theta = x_cl
	
	p = a * (1 - e^2)
	f = e * cos(omega + Omega)
	g = e * sin(omega + Omega)
	h = tan(i/2) * cos(Omega)
	k = tan(i/2) * sin(Omega)
	L = Omega + omega + theta


	return [p, f, g, h, k, L]
end

function classical_to_state_vector(x_cl, dp::DynamicsParameters)
	"""
	Reference: Fundamentals of Spacecraft Attitude Determination and Control, pg 380
	"""
	a, e, i, omega, Omega, theta = x_cl

	A11 = cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i)
	A12 = sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i)
	A13 = sin(omega) * sin(i)

	A21 = -cos(Omega) * sin(omega) - sin(Omega) * cos(omega) * cos(i)
	A22 = -sin(Omega) * sin(omega) + cos(Omega) * cos(omega) * cos(i)
	A23 = cos(omega) * sin(i)

	A = [A11 A21; A12 A22; A13 A23]

	n = sqrt(dp.mu/a^3)
	E = E_from_theta(theta, e)
	
	r_mag = a*(1 - e*cos(E))
	x_peri = a*(cos(E) - e)
	y_peri = a * sqrt(1-e^2) * sin(E)
	x_dot_peri = -(n*a^2/r_mag) * sin(E)
	y_dot_peri = (n*a^2/r_mag) * sqrt(1-e^2) * cos(E)

	r_vec = A * [x_peri; y_peri]

	v_vec = A * [x_dot_peri; y_dot_peri]

	return vcat(r_vec, v_vec)
end

function equinoctial_to_state_vector(x_eq, dp::DynamicsParameters)
	"""
	Reference: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
	"""
	p, f, g, h, k, L = x_eq

	alpha_sq = h^2 - k^2
	s_sq = 1 + h^2 + k^2
	w = 1 + f*cos(L) + g*sin(L)
	r = p/w

	r_x = (r/s_sq) * (cos(L) + alpha_sq * cos(L) + 2*h*k*sin(L))
	r_y = (r/s_sq) * (sin(L) - alpha_sq * sin(L) + 2*h*k*cos(L))
	r_z = (2*r/s_sq) * (h*sin(L) - k*cos(L))

	v_x = (-1/s_sq) * sqrt(dp.mu/p) * (sin(L) + alpha_sq*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + alpha_sq*g)
	v_y = (-1/s_sq) * sqrt(dp.mu/p) * (-cos(L) + alpha_sq*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + alpha_sq*f)
	v_z = (2/s_sq) * sqrt(dp.mu/p) * (h*cos(L) + k*sin(L) + f*h + g*k)

	return [r_x, r_y, r_z, v_x, v_y, v_z]
end

function state_vector_to_equinoctial_elements(x_st, dp::DynamicsParameters)
	"""
	Reference: Fundamentals of Spacecraft Attitude Determination and Control, pg 380
	"""

	# need better method for this
	x_cl = state_vector_to_classical_elements(x_st, dp)
	x_eq = classical_to_equinoctial_elements(x_cl)

	return x_eq
end

function state_vector_to_classical_elements(x_st, dp::DynamicsParameters)
	"""
	Source: Followed https://github.com/sisl/SatelliteDynamics.jl/blob/master/src/astrodynamics.jl#L248
	"""

	r = x_st[1:3]
	v = x_st[4:6]

	h = cross(r, v)
	W = h/norm(h)

	i = atan(sqrt(W[1]^2 + W[2]^2), W[3])
	Omega = atan(W[1], -W[2])
	p = norm(h)^2/dp.mu
	a = 1.0/(2.0/norm(r) - norm(v)^2/dp.mu)

	# numerical stability hack for circular/near circular orbits
	# ensures that (1-p/a) is always positive
	if isapprox(a, p, atol=1e-9, rtol=1e-8)
		p = a
	end

	n = sqrt(dp.mu/(a^3)) # mean motion
	e = sqrt(1 - p/a) # eccentricity
	E = atan(dot(r, v)/(n*a^2), (1 - norm(r)/a)) # eccentric anomaly
	u = atan(r[3], -r[1]*W[2] + r[2]*W[1]) # mean longitude
	theta = atan(sqrt(1-e^2)*sin(E), cos(E) - e) # true anomaly
	omega = u - theta # argument of perigee

	omega = mod2pi(omega)
	Omega = mod2pi(Omega)
	theta = mod2pi(theta)

	return [a, e, i, omega, Omega, theta]
end

function E_from_M(M, e; tol=1e-10)
    f(E) = E - e * sin(E) - M # = 0
    f_prime(E) = 1 - e * cos(E)
    
    # solve for E using Newton's method
    E_n = M
    while abs( f(E_n) ) > tol
        E_n = E_n - f(E_n)/f_prime(E_n)
    end
    
    return E_n
end

function theta_from_E(E, e)
    return atan(sqrt(1-e^2) * sin(E), cos(E) - e)
end

function E_from_theta(theta, e)
    return atan(sqrt(1 - e^2) * sin(theta), cos(theta) + e)
end

function M_from_E(E, e)
    return E - e * sin(E)
end
