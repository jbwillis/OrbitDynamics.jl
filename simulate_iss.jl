# Test/Example program of using this package to simulate the ISS
using Pkg;Pkg.activate(".");Pkg.instantiate()

include("orbit_dynamics.jl")
include("orbit_representations.jl")
include("orbit_plotting.jl")

dp = DynamicsParameters(m_satellite=1.0, A=0.1)

sma_iss = 420e3 + dp.R_earth
# e_iss = 0.0003836
e_iss = 0.1
i_iss = deg2rad(51.64)
ω_iss = deg2rad(90)
Ω_iss = deg2rad(-1)
θ0_iss = deg2rad(0.0)

x0_cl = [sma_iss, e_iss, i_iss, ω_iss, Ω_iss, θ0_iss]
x0_st = classical_to_state_vector(x0_cl, dp)
x0_eq = classical_to_equinoctial_elements(x0_cl)

t_sim = 12*Int(hr)
# solve all with 1sec timestep
x_st, t_st = solve_orbit_dynamics_ECI_state(x0_st, dp, t_sim, 1)
x_cl, t_cl = solve_orbit_dynamics_classical_elements(x0_cl, dp, t_sim, 1)
x_eq, t_eq = solve_orbit_dynamics_equinoctial(x0_eq, dp, t_sim, 1)

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
plot_classical_elements(x_cl, t_cl; fig=fig, ax=ax, label="Classical")
plot_classical_elements(cl_eq, t_eq; fig=fig, ax=ax, label="Equinoctial")

# plot specific mechanical energy
sme_st = reduce(hcat, [specific_mechanical_energy(x_st[:,i], cl_st[:,i], dp) for i=1:size(x_st)[2]])
sme_cl = reduce(hcat, [specific_mechanical_energy(rv_cl[:,i], x_cl[:,i], dp) for i=1:size(rv_cl)[2]])
sme_eq = reduce(hcat, [specific_mechanical_energy(rv_eq[:,i], cl_eq[:,i], dp) for i=1:size(rv_eq)[2]])

fig, ax = plt.subplots(1)
ax.plot(t_st ./ Int(hr), sme_st', label="Cartesian")
ax.plot(t_cl ./ Int(hr), sme_cl', label="Classical")
ax.plot(t_eq ./ Int(hr), sme_eq', "--", label="Equinoctial")
ax.set_title("Specific Mechanical Energy")
ax.legend()

# plot gravity perturbation norm
g_st = reduce(hcat, [norm(gravity_perturbation_ECI(rv_eq[:,i], dp)) for i=1:size(rv_eq)[2]])
g_cl = reduce(hcat, [norm(gravity_perturbation_classical(cl_eq[:,i], dp)) for i=1:size(cl_eq)[2]])
g_eq = reduce(hcat, [norm(gravity_perturbation_equinoctial(x_eq[:,i], dp)) for i=1:size(x_eq)[2]])

fig, ax = plt.subplots(1)
ax.plot(t_eq ./ Int(hr), g_st', label="Cartesian")
ax.plot(t_eq ./ Int(hr), g_cl', label="Classical")
ax.plot(t_eq ./ Int(hr), g_eq', "--", label="Equinoctial")
ax.legend()
ax.set_title("J2 Gravity Norm - same state, different methods")

# plot position and velocity error
fig, ax = plt.subplots(2, 1)
r_err_cl = [ norm(rv_cl[1:3,i] - x_st[1:3,i]) for i=1:size(x_st,2)]
r_err_eq = [ norm(rv_eq[1:3,i] - x_st[1:3,i]) for i=1:size(x_st,2)]

v_err_cl = [ norm(rv_cl[4:6,i] - x_st[4:6,i]) for i=1:size(x_st,2)]
v_err_eq = [ norm(rv_eq[4:6,i] - x_st[4:6,i]) for i=1:size(x_st,2)]

ax[1].plot(t_cl ./ Int(hr), r_err_cl, label="Classical")
ax[1].plot(t_eq ./ Int(hr), r_err_eq, label="Equinoctial")
ax[1].set_ylabel("Position Error")
ax[1].legend()

ax[2].plot(t_cl ./ Int(hr), v_err_cl, label="Classical")
ax[2].plot(t_eq ./ Int(hr), v_err_eq, label="Equinoctial")
ax[2].set_ylabel("Velocity Error")
ax[2].set_xlabel("Time (hr)")
fig.suptitle("Position and Velocity Errors")

