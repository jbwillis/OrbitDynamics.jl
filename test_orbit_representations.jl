using Test

include("orbit_dynamics.jl")
include("orbit_representations.jl")

@testset begin
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


