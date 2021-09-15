# include("preamble.jl")

import Pkg
Pkg.activate("../../")
Pkg.status()
Pkg.instantiate()

using Plots
using BenchmarkTools

include("csvs.jl")
include("odes.jl")
include("objectives.jl")

include("../../src/losses.jl")
include("../../src/models/ina.jl");
include("../../src/models/currents.jl")

##

kwargs_loss_robust = (α=-2, c=0.1)

p_dict = Dict{String, Float64}(zip(legend_constants.name, legend_constants.value));
p_dict["α"] = kwargs_loss_robust.α
p_dict["c"] = kwargs_loss_robust.c


## p_keys_opt

p_keys_opt_m = ["a0_m", "s_m", "b0_m", "delta_m", "v_half_m", "k_m"]
p_keys_opt_h = ["a0_h", "s_h", "b0_h", "delta_h", "v_half_h", "k_h"]
p_keys_opt_j = ["tau_j_const", "a0_j", "s_j", "b0_j", "delta_j"]

p_keys_opt_INa = ["g_max"]

p_keys_opt_other = ["α"]  # loss_robust
p_keys_opt_patch = ["c_p","c_m","R","R_f","g_leak","tau_z","x_c_comp","x_r_comp"]
p_keys_alpha = ["alpha"]

p_keys_opt = vcat(p_keys_opt_INa, p_keys_opt_other, p_keys_opt_m, p_keys_opt_h, p_keys_opt_j)

legend_subset = legend_constants[in(p_keys_opt).(legend_constants.name), :]
legend_subset.is_log .= true

mask_log = Dict(zip(legend_subset.name,
                    Array{Bool}(legend_subset.is_log)));
mask_log["α"] = Bool(0)
mask_log["c"] = Bool(1)


bounds = []

for key ∈ p_keys_opt
    if key == "α"
        append!(bounds, [(-5., 0.)])
        continue
    end
    idx = findfirst(legend_subset.name .== key)
    lb, ub = legend_subset.bound_1[idx], legend_subset.bound_2[idx]
    b = (lb, ub)
    if mask_log[key]
        b = map(x -> log.(x), b)
    end
    append!(bounds, [b])
end

bounds = Vector{Tuple{Float64, Float64}}(bounds)


## 

x₀ = [mask_log[k] ? log.(p_dict[k]) : p_dict[k] for k in p_keys_opt]
p_kwargs = (; p_keys_opt, p_dict, mask_log, kwargs_loss_robust)
p = prepare_p(x₀, p_kwargs...)

# rhs = ODEFunction(compute_rates!, syms=[:v_comp, :v_p, :v_m, :m, :h, :j, :I_out])
rhs = ODEFunction((du, u, p, t) -> compute_rates!(du, u, p, t; 1e-6),
                  syms=[:v_comp, :v_p, :v_m, :m, :h, :j, :I_out]);

cb_step_v1  = PresetTimeCallback(protocol.t, change_step_v1!, save_positions=(false, false))
prob = ODEProblem(rhs, u₀, tspan, p, callback=cb_step_v1)
solve_kwargs_default = (; reltol, abstol, solver, saveat, dt)


##

sol = solve_model(prob, (;), solve_kwargs_default);
data_true = similar(sol[:I_out])
data_buffer = similar(data_true)

calculate_data_segments!(data_true, p, prob, solve_kwargs_default)
calculate_loss_segments!(data_buffer, data_true, x₀, p_kwargs, prob, solve_kwargs_default)
