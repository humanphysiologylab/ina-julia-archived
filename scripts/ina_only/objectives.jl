using Parameters


function prepare_p(x, p_keys_opt, p_dict, mask_log, kwargs_loss)
    
    @assert length(x) == length(p_keys_opt)
    
    p = deepcopy(p_dict)
    for (k, v, is_log) in zip(p_keys_opt, x, mask_log)
        p[k] = mask_log[k] ? exp(v) : v
    end
    
    p["α"] = kwargs_loss.α
    
    return p
    
end


function sanity_checks(p)  # TODO
    v_m_range = -100:1:50
    
    τ_m, τ_h, τ_j = [map(v_m -> calculate_τ(v_m, p), v_m_range) for calculate_τ ∈ (calculate_tau_m,
                                                                                   calculate_tau_h, 
                                                                                   calculate_tau_j)]
    
    (τ_m_max, τ_m_idxmax), (τ_h_max, τ_h_idxmax), (τ_j_max, τ_j_idxmax) = map(τ -> findmax(τ), (τ_m, τ_h, τ_j))
    τ_m_min, τ_h_min, τ_j_min = map(τ -> min(τ...), (τ_m, τ_h, τ_j))

    flags = [all(1e-7 .< τ_m_min .< 1e-3),
             all(1e-8 .< τ_h_min .< 1e-2),
             all(5e-4 .< τ_j_min .< 2.),
             τ_m_max > 4e-5,
             τ_h_max > 1e-3,
             τ_j_max > 8e-2,
             (-90. < v_m_range[τ_m_idxmax] < 0.),
             (-90. < v_m_range[τ_h_idxmax] < 0.),
             (-90. < v_m_range[τ_j_idxmax] < 0.)]
end


function calculate_u∞!(u∞, p, prob, solve_kwargs_default)::Symbol  # returns status

    @unpack reltol, abstol, solver, saveat, dt = solve_kwargs_default

    t_eq_end = 10.
    remake_kwargs = (p=p, tspan=(0., t_eq_end), callback=nothing)
     solve_kwargs = (saveat=[t_eq_end], reltol, abstol, solver, dt)

    sol_eq = solve_model(prob, remake_kwargs, solve_kwargs)
    u∞[:] = sol_eq.u[1]
    return sol_eq.retcode

end


function calculate_data_segments_parallel(data_segmented, u∞, p, prob, solve_kwargs_default, n_steps=20)::Symbol  # returns status

    tspan_end = 5.
    step_size = 5e-5
    segment_size_sec = tspan_end / n_steps
    segment_size_timestamps = Int(tspan_end / n_steps / step_size)

    Threads.@threads for i ∈ 1: n_steps

        tspan_segment = (segment_size_sec * (i - 1),
                         segment_size_sec * i)
        remake_kwargs = (u0=u∞, p=deepcopy(p), tspan=tspan_segment)
        sol_segment = solve_model(prob, remake_kwargs, solve_kwargs_default)

        if sol_segment.retcode == :Success
            i_start = (i - 1) * segment_size_timestamps + 1
            i_end = i * segment_size_timestamps + 1
            data_segmented[i_start: i_end] = sol_segment[:I_out]
        else
            sol_segment.retcode
        end

    end

    return :Success

end


function calculate_data_segments_serial(data_segmented, p, prob, solve_kwargs_default, n_steps=20)::Symbol  # returns status

    tspan_end = 5.
    step_size = 5e-5
    segment_size_sec = tspan_end / n_steps
    segment_size_timestamps = Int(tspan_end / n_steps / step_size)

    for i ∈ 1: n_steps

        tspan_segment = (segment_size_sec * (i - 1),
                         segment_size_sec * i)
        remake_kwargs = (u0=u∞, p=deepcopy(p), tspan=tspan_segment)
        sol_segment = solve_model(prob, remake_kwargs, solve_kwargs_default)

        if sol_segment.retcode == :Success
            i_start = (i - 1) * segment_size_timestamps + 1
            i_end = i * segment_size_timestamps + 1
            data_segmented[i_start: i_end] = sol_segment[:I_out]
        else
            sol_segment.retcode
        end

    end

    return :Success

end


function calculate_data_segmented!(data_segmented, p, prob, solve_kwargs_default, n_steps=20)::Symbol

    u∞ = similar(prob.u0)
    retcode = calculate_u∞!(u∞, p, prob, solve_kwargs_default)

    if retcode != :Success
        return retcode
    end

    retcode = calculate_data_segments_parallel(data_segmented, u∞, p, prob, solve_kwargs_default, n_steps)

    if retcode != :Success
        return retcode
    end

    return :Success

end


function calculate_loss_segments!(data_segmented, data_true,
                                 x, p_kwargs, prob, solve_kwargs_default, n_steps=20)::Float64

    p = prepare_p(x, p_kwargs...)
    is_ok = all(sanity_checks(p))
    if !is_ok
        return Inf
    end    

    retcode = calculate_data_segmented!(data_segmented, p, prob, solve_kwargs_default, n_steps)

    if retcode == :Success
        residuals = data_segmented - data_true
        @unpack c, α = p
        loss = calculate_loss_robust.(residuals, α=α, c=c)
        loss = mean(loss)
    else
        loss = Inf
    end
    
    return loss
        
end