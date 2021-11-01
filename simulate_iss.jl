# Test/Example program of using this package to simulate the ISS

include("orbit_dynamics.jl")
include("orbit_representations.jl")
include("orbit_plotting.jl")

dp = DynamicsParameters(m_satellite=1.0, A = .1, J2=0.0, rho=0.0)
sma_iss = 420e3 + dp.R_earth
i_iss = deg2rad(51.64)
e_iss = 0.1
ω_iss = deg2rad(90)
Ω_iss = deg2rad(0)
θ0_iss = deg2rad(0)

x0_cl = [sma_iss, i_iss, e_iss, ω_iss, Ω_iss, θ0_iss]
x0_st = classical_to_state_vector(x0_cl, dp)
x0_eq = classical_to_equinoctial_elements(x0_cl)

t_sim = 12*Int(hr)
x_st, t_st = solve_orbit_dynamics_ECI_state(x0_st, dp, t_sim)
x_cl, t_cl = solve_orbit_dynamics_classical_elements(x0_cl, dp, t_sim)
x_eq, t_eq = solve_orbit_dynamics_equinoctial(x0_eq, dp, t_sim)

# plot position in 3D
rv_cl = reduce(hcat, [classical_to_state_vector(x_cl[:,i], dp) for i=1:size(x_cl)[2]])
rv_eq = reduce(hcat, [equinoctial_to_state_vector(x_eq[:,i], dp) for i=1:size(x_eq)[2]])

fig_3d, ax_3d = plot_3D_position(x_st[1:3,:]; label="Cartesian")
plot_3D_position(rv_cl[1:3,:]; fig=fig_3d, ax=ax_3d, label="Classical")
plot_3D_position(rv_eq[1:3,:]; fig=fig_3d, ax=ax_3d, label="Equinoctial")

# plot orbital elements
cl_st = reduce(hcat, [state_vector_to_classical_elements(x_st[:,i], dp) for i=1:size(x_st)[2]])
cl_eq = reduce(hcat, [equinoctial_to_classical_elements(x_eq[:,i]) for i=1:size(x_eq)[2]])

fig, ax = plot_classical_elements(cl_st, t_st; label="Cartesian")
plot_classical_elements(x_cl, t_eq; fig=fig, ax=ax, label="Classical")
plot_classical_elements(cl_eq, t_eq; fig=fig, ax=ax, label="Equinoctial")








