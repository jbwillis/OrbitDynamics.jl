@testset "Gravity, Drag in Classical and Equinoctial" begin

	dp = DynamicsParameters(m_satellite_si=1.0, A_si=10)
	sma_iss = 420e3 + dp.R_earth
	e_iss = 0.1
	i_iss = deg2rad(51.64)
	ω_iss = deg2rad(90)
	Ω_iss = deg2rad(0)
	θ0_iss = deg2rad(45)

	x0_cl = [sma_iss, e_iss, i_iss, ω_iss, Ω_iss, θ0_iss]
	x0_eq = classical_to_equinoctial_elements(x0_cl)

	d_eq = drag_perturbation_equinoctial(x0_eq, dp)
	d_cl = drag_perturbation_classical(x0_cl, dp)

	@test d_eq ≈ d_cl

	g_eq = gravity_perturbation_equinoctial(x0_eq, dp)
	g_cl = gravity_perturbation_classical(x0_cl, dp)

	@test g_eq ≈ g_cl

end

