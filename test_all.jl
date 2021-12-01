using Test

include("orbit_dynamics.jl")
include("orbit_representations.jl")

@testset "Representation conversions are identity" begin
	dp = DynamicsParameters()
	x_cl_1 = [1e4 + dp.R_earth, # a
			0.5, # e
			deg2rad(45), # i 
			deg2rad(90), # omega
			deg2rad(125), # Omega
			deg2rad(15)] # theta
	x_eq_1 = classical_to_equinoctial_elements(x_cl_1)
	x_st_1 = classical_to_state_vector(x_cl_1, dp)

	@test x_cl_1 ≈ equinoctial_to_classical_elements(x_eq_1)
	@test x_cl_1 ≈ state_vector_to_classical_elements(x_st_1, dp)
	@test x_eq_1 ≈ state_vector_to_equinoctial_elements(x_st_1, dp)
	@test x_st_1 ≈ equinoctial_to_state_vector(x_eq_1, dp)

end

@testset "M, E, theta conversions are identity" begin
	# M -> E -> theta -> E -> M
	e = 0.3
    M0 = pi/3
    E1 = E_from_M(M0, e)
    M1 = M_from_E(E1, e)
    theta = theta_from_E(E1, e)
    E2 = E_from_theta(theta, e)
    M2 = M_from_E(E1, e)
    
    @test abs(M0 - M1) < 1e-6
    @test abs(E1 - E2) < 1e-6
    @test abs(M0 - M2) < 1e-6
end

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

@testset "Scaling produces similar results" begin

	dp_scaled = DynamicsParameters(distance_scale=1e7, time_scale=3600)
	dp_unscaled = DynamicsParameters(distance_scale=1.0, time_scale=1.0)

	x0_unscaled = classical_to_state_vector([dp_unscaled.R_earth + 500e3, 0.1, deg2rad(45), deg2rad(60), deg2rad(240), deg2rad(145)], dp_unscaled)
	x0_scaled = scale_state_vector(x0_unscaled, dp_scaled)
	@test x0_unscaled ≈ unscale_state_vector(x0_scaled, dp_scaled)
	@test x0_scaled ≈ scale_state_vector(unscale_state_vector(x0_scaled, dp_scaled), dp_scaled)

	x0_dot_unscaled = orbit_dynamics_ECI_state(x0_unscaled, dp_unscaled, 0.0)
	x0_dot_scaled = orbit_dynamics_ECI_state(x0_scaled, dp_scaled, 0.0)

	@test x0_dot_scaled ≈ scale_state_vector_dot(x0_dot_unscaled, dp_scaled)
	@test x0_dot_unscaled ≈ unscale_state_vector_dot(x0_dot_scaled, dp_scaled)

end

