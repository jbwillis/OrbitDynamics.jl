using PyPlot
using LaTeXStrings

@enum t_scale_e sec=1 min=60 hr=60*60 day=24*60*60 year=365*24*60*60


function plot_3D_position(r; label=nothing, fig=nothing, ax=nothing)
	if fig == nothing || ax == nothing
		fig = plt.figure()

		# 3D in top right
		ax3D = fig.add_subplot(2, 2, 2, projection="3d")
		ax3D.scatter3D(0.0, 0.0, 0.0, color="green")

		# 2D in other quadrants
		ax_ij = fig.add_subplot(2, 2, 1)
		ax_ij.scatter(0.0, 0.0, color="green")
		ax_ij.set_xlabel("I")
		ax_ij.set_ylabel("J")

		ax_ik = fig.add_subplot(2, 2, 3)
		ax_ik.scatter(0.0, 0.0, color="green")
		ax_ik.set_xlabel("I")
		ax_ik.set_ylabel("K")

		ax_jk = fig.add_subplot(2, 2, 4)
		ax_jk.scatter(0.0, 0.0, color="green")
		ax_jk.set_xlabel("J")
		ax_jk.set_ylabel("K")

		ax = [ax_ij, ax3D, ax_ik, ax_jk]
	end
	ax[2].plot3D(r[1,:], r[2,:], r[3,:], label=label)
    ax[2].set_title("3D Position")
    ax[2].set_xlabel("I")
    ax[2].set_ylabel("J")
    ax[2].set_zlabel("K")

	ax[1].plot(r[1,:], r[2,:])
	ax[3].plot(r[1,:], r[3,:])
	ax[4].plot(r[2,:], r[3,:])

	ax[1].axis("equal")
	ax[3].axis("equal")
	ax[4].axis("equal")

	if label != nothing
		ax[2].legend()
	end

	return fig, ax
end

function plot_classical_elements(x_cl, t; label=nothing, fig=nothing, ax=nothing, t_scale::t_scale_e=hr)
	if fig == nothing || ax == nothing
		fig, ax = plt.subplots(6, 1, sharex=true)

		fig.suptitle("Osculating Elements")
	end

	t_plot = t ./ Int(t_scale)
	ax[1].plot(t_plot, x_cl[1,:] ./ 1e3, label=label)
	ax[1].set_ylabel(L"SMA (km)")

	ax[2].plot(t_plot, x_cl[2,:])
	ax[2].set_ylabel(L"$e$")

	ax[3].plot(t_plot, rad2deg.(mod2pi.(x_cl[3,:])))
	ax[3].set_ylabel(L"$i$ (deg)")

	ax[4].plot(t_plot, rad2deg.(mod2pi.(x_cl[4,:])))
	ax[4].set_ylabel(L"$\omega$ (deg)")

	ax[5].plot(t_plot, rad2deg.(mod2pi.(x_cl[5,:])))
	ax[5].set_ylabel(L"$\Omega$ (deg)")

	ax[6].plot(t_plot, rad2deg.(mod2pi.(x_cl[6,:])))
	ax[6].set_ylabel(L"$\theta$ (deg)")

	ax[6].set_xlabel(L"Time (" * string(t_scale) * L")")

	if label != nothing
		ax[1].legend()
	end

	return fig, ax
end


