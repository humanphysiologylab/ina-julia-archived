using Parameters


function calculate_v_cp(v_comp, p)
    @unpack v_c, alpha = p
    v_cp = v_c + (v_c - v_comp) * (1 / (1 - alpha) - 1)
end


function calculate_d_v_comp(v_comp, p)
    @unpack v_c, x_c_comp, c_m, x_r_comp, R, alpha = p
    d_v_comp = (v_c - v_comp) / (x_c_comp * c_m * x_r_comp * R * (1 - alpha))
end


function calculate_d_v_p(v_cp, v_p, p)
    @unpack c_p, R_f = p
    d_v_p = (v_cp - v_p) / (c_p * R_f)
end


function calculate_d_v_m(v_m, v_p, I_leak, I_Na, p)
    @unpack v_off, R, c_m = p
    d_v_m = (v_p + v_off - v_m) / (R * c_m) - 1e-9 * (I_leak + I_Na) / c_m
end


function calculate_I_leak(v_m, p)
    @unpack g_leak = p
    I_leak = g_leak * v_m
end


function calculate_I_c(v_m, v_p, I_leak, I_Na, p)
    @unpack c_m = p
    I_c = 1e9 * p.c_m * calculate_d_v_m(v_m, v_p, I_leak, I_Na, p)
end


function calculate_I_p(v_cp, v_p, p)
    @unpack c_p = p
    I_p = 1e9 * c_p * calculate_d_v_p(v_cp, v_p, p)
end


function calculate_I_comp(v_comp, p)
    @unpack v_c, x_c_comp, c_m = p
    I_comp = 1e9 * x_c_comp * c_m * calculate_d_v_comp(v_comp, p)
end


function calculate_d_I_out(I_in, I_out, p)
    @unpack tau_z = p
    d_I_out = (I_in - I_out) / tau_z
end
