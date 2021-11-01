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

	dp = DynamicsParameters(m_satellite=1.0, A=10)
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
