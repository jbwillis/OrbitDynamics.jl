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

@testset "Scaling produces similar results" begin

	dp_scaled = DynamicsParameters(distance_scale=1e7, time_scale=3600)
	dp_unscaled = DynamicsParameters(distance_scale=1.0, time_scale=1.0)

	x0_unscaled = classical_to_state_vector([dp_unscaled.R_earth + 500e3, 0.1, deg2rad(45), deg2rad(60), deg2rad(240), deg2rad(145)], dp_unscaled)
	x0_scaled = scale_state_vector(x0_unscaled, dp_scaled)
	@test x0_unscaled ≈ unscale_state_vector(x0_scaled, dp_scaled)
	@test x0_scaled ≈ scale_state_vector(unscale_state_vector(x0_scaled, dp_scaled), dp_scaled)

    x0_cl_1 = state_vector_to_classical_elements(x0_scaled, dp_scaled)
    x0_cl_2 = state_vector_to_classical_elements(x0_unscaled, dp_unscaled)
	@test x0_cl_1[2:end] ≈ x0_cl_2[2:end] # scaling shouldn't affect angles
	@test 0 < x0_cl_1[1] < 1.0
	@test 1e6 < x0_cl_2[1] < 1e7

	x0_cl_1_cart = classical_to_state_vector(x0_cl_1, dp_scaled)
	@test x0_cl_1_cart ≈ x0_scaled
	x0_cl_2_cart = classical_to_state_vector(x0_cl_2, dp_unscaled)
	@test x0_cl_2_cart ≈ x0_unscaled

	x0_dot_unscaled = orbit_dynamics_ECI_state(x0_unscaled, dp_unscaled, 0.0)
	x0_dot_scaled = orbit_dynamics_ECI_state(x0_scaled, dp_scaled, 0.0)

	@test x0_dot_scaled ≈ scale_state_vector_dot(x0_dot_unscaled, dp_scaled)
	@test x0_dot_unscaled ≈ unscale_state_vector_dot(x0_dot_scaled, dp_scaled)

end

@testset "Cartesian to cylindrical" begin
@test pos_cartesian_to_cylindrical(pos_cylindrical_to_cartesian([1,2,3])) ≈ [1, 2, 3]
@test vel_cartesian_to_cylindrical(pos_cylindrical_to_cartesian([5, 1, -1]), vel_cylindrical_to_cartesian([5, 1, -1], [3, 0.2, 1])) ≈ [3, 0.2, 1]
@test accel_cylindrical_to_cartesian(pos_cartesian_to_cylindrical([1, 0, 0]), vel_cartesian_to_cylindrical([1, 0, 0], [9, 8, 7]), accel_cartesian_to_cylindrical([1, 0, 0], [9, 8, 7], [11, 12, 13])) ≈ [11, 12, 13]

dp = DynamicsParameters()
x0_cl = [420e3 + dp.R_earth, 0.1, deg2rad(45), deg2rad(155), deg2rad(247), deg2rad(101)]
xu0 = [classical_to_state_vector(x0_cl, dp); 1.0]
xu0_cyl = [state_cartesian_to_cylindrical(xu0[1:6]); xu0[7]]

x_cart_dot = orbit_dynamics_ECI_drag_control(xu0, dp, 0.0)
x_cyl_dot = orbit_dynamics_cylindrical_drag_control(xu0_cyl, dp, 0.0)

# make sure dynamics are the same in cartesian and cylindrical
@test x_cart_dot[1:6] ≈ state_dot_cylindrical_to_cartesian(xu0_cyl, x_cyl_dot)
@test x_cyl_dot[1:6] ≈ state_dot_cartesian_to_cylindrical(xu0, x_cart_dot)

end
