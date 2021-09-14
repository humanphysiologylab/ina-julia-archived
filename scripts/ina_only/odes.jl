using DifferentialEquations
using Sundials


const u₀ = [-80., -80., -80., 0., 1., 1., 0.];
const reltol = 1e-3
const abstol = [1e-2, 1e-2, 1e-2, 1e-4, 1e-4, 1e-4, 1e-2]
const solver = CVODE_BDF();
const dt = 1e-9  # initial
const tspan_initial = (0., 10.)
const tspan = (0., 5.)
const saveat = 5e-5 # tspan[1]: 5e-5: tspan[2]

include("../../src/models/ina.jl");
include("../../src/models/currents.jl");

function compute_algebraic(du, u, p, t)
    
    v_comp, v_p, v_m, m, h, j, I_out = u
        
    tau_m  = calculate_tau_m(v_m, p)
    tau_h  = calculate_tau_h(v_m, p)
    tau_j  = calculate_tau_j(v_m, p)
    
    m_inf  = calculate_m_inf(v_m, p)
    h_inf  = calculate_h_inf(v_m, p)
    
    v_cp   = calculate_v_cp(v_comp, p)

    I_leak = 0 # calculate_I_leak(v_m, p)
    I_Na   = calculate_I_Na(v_m, m, h, j, p)
    I_c    = 0 # calculate_I_c(v_m, v_cp, I_leak, I_Na,  p)  # or calculate_I_c(v_m, v_p, I_leak, I_Na,  p)
    I_p    = 0  # or calculate_I_p(v_cp, v_p, p)
    I_comp = 0 # calculate_I_comp(v_comp, p)
    
    I_in   = I_leak + I_Na + I_c - I_comp + I_p
    
    a = (;
     tau_m, tau_h, tau_j, m_inf, h_inf,
     v_cp,
     I_leak, I_Na, I_c, I_comp, I_in
    )
    
end


function compute_rates!(du, u, p, t; safety_factor=0.)
    
    v_comp, v_p, v_m, m, h, j, I_out = u
    
    a = compute_algebraic(du, u, p, t)
    
    du[1] = calculate_d_v_comp(v_comp, p)  # v_comp
    
    @unpack v_cp, I_leak, I_Na = a
    du[2] = 0#calculate_d_v_p(v_cp, v_p, p)  # v_p
    du[3] = calculate_d_v_m(v_m, v_cp, I_leak, I_Na, p)  # v_m
        
    @unpack m_inf, tau_m, h_inf, tau_h, tau_j = a
    du[4] = calculate_d_gate(m_inf, m, tau_m + safety_factor)  # m
    du[5] = calculate_d_gate(h_inf, h, tau_h + safety_factor)  # h
    du[6] = calculate_d_gate(h_inf, j, tau_j)  # j
        
    @unpack I_in = a
    @unpack tau_z = p
    du[7] = calculate_d_gate(I_in, I_out, tau_z)  # I_out
    
    nothing
end


function change_step_v1!(integrator)
    t = integrator.t
    v_c = find_step(t)
    integrator.p["v_c"] = v_c
    set_proposed_dt!(integrator, 1e-7)
    nothing
end


function solve_model(prob, remake_kwargs, solve_kwargs)
    prob_remade = remake(prob; remake_kwargs...)
    sol = solve(prob_remade; solve_kwargs...)
end
