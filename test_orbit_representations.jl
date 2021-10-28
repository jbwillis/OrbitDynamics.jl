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


