module SingleOrbitTools
using LinearAlgebra
using OrbitDynamics

export one_orbit
"""
    x_T = one_orbit(x0; dynamics=orbit_dynamics_ECI_state, plane_normal=nothing, disp=false, tol=1e-10, return_all=false, return_T=false, dp=DynamicsParameters(), t_step=1.0)

Integrate orbit dynamics over one orbit.
Can use different dynamics, but first 6 elements must be ECI position and velocity
"""
function one_orbit(x0; 
        dynamics=orbit_dynamics_ECI_state, 
        plane_normal=nothing, 
        disp=false, 
        tol=1e-10, 
        return_all=false, 
        return_T=false, 
        dp=DynamicsParameters(), 
        t_step=1.0)
    t = 0.0

    if isnothing(plane_normal)
        # use tangent vector at x0 as vector defining transverse plane normal
        r0 = x0[1:3]
        v0 = x0[4:6]
        h = cross(r0, v0)
        T0 = cross(h, r0)
        T0 /= norm(T0)
        plane_normal = T0
    end
    
    
    x0_cl = state_vector_to_classical_elements(x0[1:6], dp)
    orbital_period = 2*pi * sqrt(x0_cl[1]^3 / dp.mu)

    # initial conditions
    x_km1 = x0
    x_k = x0
    crossed_plane = false
    
    if return_all
        x_all = []
        t_all = []
        push!(x_all, x0)
        push!(t_all, t)
    end

    # integrate until crossed the plane defined by plane_normal
    while !crossed_plane

        # update
        x_km1 = x_k
        x_k = step_RK4(dynamics, x_k, dp, t, t_step)
        t += t_step
        
        if return_all
            push!(x_all, x_k)
            push!(t_all, t)
        end
        
        crossed_plane = (t > orbital_period/2) && (x_km1[1:3]' * plane_normal < 0) && (x_k[1:3]' * plane_normal >= 0)
    end
    
    # now perform a binary search to refine 
    t_search = - t_step/2
    bi = 0
    while abs(x_k[1:3]' * plane_normal) > tol
        # update
        x_km1 = x_k
        x_k = step_RK4(dynamics, x_k, dp, t, t_search)
        t += t_search
        
        if (x_k[1:3]' * plane_normal >= 0)
            # still on same side of plane - keep going backwards
            t_search = -abs(t_search)/2
        else
            # jumped across plane - go forwards
            t_search = abs(t_search)/2
        end
        
        if return_all
            push!(x_all, x_k)
            push!(t_all, t)
        end
        
        bi += 1
    end
    
    if disp
        println("Integrated time = ", t)
        println("Orbital Period = ", orbital_period)
        
        println("Binary search took $bi iterations")
        @show t_search
        println("r^T p =", x_k[1:3]' * plane_normal)
        println()
    end
    
    if return_all
        if return_T
            return x_all, t_all
        else
            return x_all
        end
    end
    
    if return_T
        return x_k, t
    end
    
    return x_k
end

export integrate_to_T
"""
Integrate dynamics to a given time T.
Uses the same binary search sequence as `one_orbit()` to ensure that the results are the same for the same `T`

Can use different dynamics, but first 6 elements must be ECI position and velocity
"""
function integrate_to_T(x0, T_end; dynamics=orbit_dynamics_ECI_state, dp=DynamicsParameters(), t_step=1.0, return_all=false)

    t = 0.0

    if return_all
        x_all = []
        t_all = []
        push!(x_all, x0)
        push!(t_all, t)
    end
    
    # initial conditions
    x_km1 = x0
    x_k = x0
    
    while t <= T_end

        # update
        x_km1 = x_k
        x_k = step_RK4(dynamics, x_k, dp, t, t_step)
        t += t_step

        if return_all
            push!(x_all, x_k)
            push!(t_all, t)
        end
    end
    
    
    t_search = - t_step/2
    while t != T_end
        # update
        x_km1 = x_k
        x_k = step_RK4(dynamics, x_k, dp, t, t_search)
        t += t_search

        if return_all
            push!(x_all, x_k)
            push!(t_all, t)
        end
        
        if T_end - t < 0.0
            # still on same side of T_end - go backwards
            t_search = -abs(t_search)/2
        else
            # jumped across T_end - go forwards
            t_search = abs(t_search)/2
        end
    end
    
    if return_all
        return x_all, t_all
    end

    return x_k, t
end

end # Module