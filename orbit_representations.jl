"""
Code for converting between different orbit representations
"""

function equinoctial_to_classical_elements(x_eq)
	p, f, g, h, k, L = x_eq

	a = p / (1 - f^2 - g^2)
	e = sqrt(f^2 + g^2)
	i = atan(2*sqrt(h^2 + k^2), 1 - h^2 - k^2)
	omega = atan(g*h - f*k, f*h + g*k)
	Omega = atan(k, h)
	theta = L - atan(g/f)

	return [a, e, i, omega, Omega, theta]
end

function classical_to_equinoctial_elements(x_cl)

	a, e, i, omega, Omega, theta = x_cl
	
	p = q * (1 - e^2)
	f = e * cos(omega + Omega)
	g = e * sin(omega + Omega)
	h = tan(i/2) * cos(Omega)
	k = tan(i/2) * sin(Omega)
	L = Omega + omega + theta


	return [p, f, g, h, k, L]

end

function classical_to_state_vector(x_cl)

	a, e, i, omega, Omega, theta = x_cl
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

function state_vector_to_equinoctial_elements(x_st)

end

function state_vector_to_orbital_elements(x_st)

end

