using Parameters


function calculate_tau_m(v_m, p)
    @unpack a0_m, s_m, b0_m, delta_m = p
    tau_m = 1 / (a0_m * exp(v_m / s_m) + b0_m * exp(-v_m / delta_m))
end


function calculate_tau_h(v_m, p)
    @unpack a0_h, s_h, b0_h, delta_h = p
    tau_h = 1 / (a0_h * exp(-v_m / s_h) + b0_h * exp(v_m / delta_h))
end


function calculate_tau_j(v_m, p)
    @unpack tau_j_const, a0_j, s_j, b0_j, delta_j = p
    tau_j = tau_j_const + 1 / (a0_j * exp(-v_m / s_j) + b0_j * exp(v_m / delta_j))
end


function calculate_m_inf(v_m, p)
    @unpack v_half_m, k_m = p
    m_inf = 1 / (1 + exp(-(v_half_m + v_m) / k_m))
end


function calculate_h_inf(v_m, p)
    @unpack v_half_h, k_h = p
    h_inf = 1 / (1 + exp((v_half_h + v_m) / k_h))
end


function calculate_d_gate(gate_inf, gate, tau_gate)
    d_gate = (gate_inf - gate) / tau_gate
end


function calculate_I_Na(v_m, m, h, j, p)
    @unpack g_max, v_rev = p
    I_Na = g_max * h * m^3 * j * (v_m - v_rev)
end
