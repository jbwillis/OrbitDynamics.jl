using PyPlot
using LaTeXStrings

@enum t_scale_e sec=1 min=60 hr=60*60 day=24*60*60 year=365*24*60*60


function plot_3D_position(r; label=nothing, fig=nothing, ax=nothing)
	if fig == nothing || ax == nothing
		fig = plt.figure()
		ax = Axes3D(fig, auto_add_to_figure=false)
		ax.scatter3D(0.0, 0.0, 0.0, color="green")
	end
    ax.plot3D(r[1,:], r[2,:], r[3,:], label=label)
    ax.set_title("3D Position")
    ax.set_xlabel("I")
    ax.set_ylabel("J")
    ax.set_zlabel("K")
    fig.add_axes(ax)

	if label != nothing
		ax.legend()
	end

	return fig, ax
end

function plot_classical_elements(x_cl, t; label=nothing, fig=nothing, ax=nothing, t_scale::t_scale_e=hr)
	if fig == nothing || ax == nothing
		fig, ax = plt.subplots(6, 1, sharex=true)
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


