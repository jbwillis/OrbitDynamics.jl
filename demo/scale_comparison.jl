using OrbitDynamics
using OrbitDynamics.OrbitPlotting

dp_scaled = DynamicsParameters(distance_scale=1e7, time_scale=1e3)
dp_unscaled = DynamicsParameters(distance_scale=1.0, time_scale=1.0)

x0_unscaled = classical_to_state_vector([dp_unscaled.R_earth + 500e3, 0.1, deg2rad(90), deg2rad(0), deg2rad(0), deg2rad(0)], dp_unscaled)
x0_scaled = scale_state_vector(x0_unscaled, dp_scaled)

@show x0_unscaled
@show x0_scaled

x0_dot_unscaled = orbit_dynamics_ECI_state(x0_unscaled, dp_unscaled, 0.0)
x0_dot_scaled = orbit_dynamics_ECI_state(x0_scaled, dp_scaled, 0.0)

@show x0_dot_unscaled
@show x0_dot_scaled

x_step_unscaled = step_RK4(orbit_dynamics_ECI_state, x0_unscaled, dp_unscaled, 0.0, 1.0)
x_step_scaled = step_RK4(orbit_dynamics_ECI_state, x0_scaled, dp_scaled, 0.0, 1.0 / dp_scaled.time_scale)

@show x_step_unscaled
@show x_step_scaled

x_unscaled, t_unscaled = solve_orbit_dynamics_ECI_state(x0_unscaled, dp_unscaled, 3, 1.0)
x_scaled, t_scaled = solve_orbit_dynamics_ECI_state(x0_scaled, dp_scaled, 3 / dp_scaled.time_scale, 1.0 / dp_scaled.time_scale)

println("x_unscaled = \n", x_unscaled)
println("x_scaled = \n", x_scaled)

x_unscaled_all, t_unscaled_all = solve_orbit_dynamics_ECI_state(x0_unscaled, dp_unscaled, 5000, 1.0)
x_scaled_all, t_scaled_all = solve_orbit_dynamics_ECI_state(x0_scaled, dp_scaled, 5000 / dp_scaled.time_scale, 1.0 / dp_scaled.time_scale)


plot_3D_position(x_unscaled_all[1:3,:]; label="Unscaled")
plot_3D_position(x_scaled_all[1:3,:]; label="Scaled")
