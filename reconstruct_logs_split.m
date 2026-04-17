function log = reconstruct_logs_split(p,t,x)

log = struct();

% ============================================================
% State unpacking (vectorized)
% ============================================================% ============================================================
i1d   = x(p.state_idx.i1d,:);
i1q   = x(p.state_idx.i1q,:);
i2d   = x(p.state_idx.i2d,:);
i2q   = x(p.state_idx.i2q,:);
vcd   = x(p.state_idx.vcd,:);
vcq   = x(p.state_idx.vcq,:);
delta_g = x(p.state_idx.delta_g,:);
omega_g = x(p.state_idx.omega_g,:);
delta_conv = x(p.state_idx.delta_c,:);
omega_conv = x(p.state_idx.omega_c,:);
xi_Q  = x(p.state_idx.xi_Q,:);
xi_id = x(p.state_idx.xi_id,:);
xi_iq = x(p.state_idx.xi_iq,:);
it1d  = x(p.state_idx.it1d,:);
it1q  = x(p.state_idx.it1q,:);
vt1d  = x(p.state_idx.vt1d,:);
vt1q  = x(p.state_idx.vt1q,:);
it2d  = x(p.state_idx.it2d,:);
it2q  = x(p.state_idx.it2q,:);
vt2d  = x(p.state_idx.vt2d,:);
vt2q  = x(p.state_idx.vt2q,:);
xi_vd = x(p.state_idx.xi_vd,:);
xi_vq = x(p.state_idx.xi_vq,:);

% ============================================================
% Inputs / disturbances (time-varying)
% ============================================================
V_ref = arrayfun(p.V_ref,t);
Q_ref = arrayfun(p.Q_ref,t);
Vg_mag       = arrayfun(p.vg_mag,t);
theta_g = arrayfun(p.vg_phase_rad,t);

% absolute grid angle
delta_g_abs = delta_g + theta_g;

% ============================================================
% Grid voltage
% ============================================================
vg_d = Vg_mag .* cos(delta_g_abs);
vg_q = Vg_mag .* sin(delta_g_abs);

log.vg_d = vg_d;
log.vg_q = vg_q;

% ============================================================
% PCC voltage
% ============================================================
ilpd = i1d - i2d - it1d - it2d;
ilpq = i1q - i2q - it1q - it2q;

vPCCd = vcd + p.Rlp .* ilpd;
vPCCq = vcq + p.Rlp .* ilpq;

log.V_PCC = [vPCCd; vPCCq];
log.V_PCC_mag = sqrt(vPCCd.^2 + vPCCq.^2);

% ============================================================
% Power at PCC
% ============================================================
P_PCC = vPCCd .* i2d + vPCCq .* i2q;
Q_PCC = vPCCq .* i2d - vPCCd .* i2q;

log.P_PCC = P_PCC;
log.Q_PCC = Q_PCC;

% ============================================================
% Voltage magnitude
% ============================================================
V_meas = sqrt(vPCCd.^2 + vPCCq.^2);


% ============================================================
% Q-control (your formulation)
% ============================================================
Q_meas = vPCCq .* i2d - vPCCd .* i2q;

err_Q = Q_ref - Q_meas;
err_v = V_ref - V_meas;

Vref_cf = ...
    + p.Kpq .* err_Q ...
    + p.Kpq .* p.Kvq .* err_v ...
    + p.Kiq .* xi_Q;



log.Vrefd = Vref_cf;

% ============================================================
% dq transformation of voltage reference
% ============================================================
c = cos(delta_conv);
s = sin(delta_conv);

Vref_d = c .* Vref_cf;
Vref_q = s .* Vref_cf;

log.Vref = [Vref_d; Vref_q];

% ============================================================
% Virtual impedance dynamics (for interpretation)
% ============================================================
e_vd = Vref_d - vPCCd;
e_vq = Vref_q - vPCCq;

d_vv_d = p.Kpv .* e_vd + p.Kiv .* xi_vd;
d_vv_q = p.Kpv .* e_vq + p.Kiv .* xi_vq;

Zv2 = sqrt(p.Rv^2 + p.Xv^2);

Id_ref = (d_vv_d .* p.Rv + d_vv_q .* p.Xv) ./ Zv2;
Iq_ref = (d_vv_q .* p.Rv - d_vv_d .* p.Xv) ./ Zv2;

% ============================================================
% Current limit
% ============================================================
I_mag = sqrt(Id_ref.^2 + Iq_ref.^2);

scale = ones(size(I_mag));

idx = I_mag > p.I_max;

scale(idx) = p.I_max ./ I_mag(idx);

Id_ref_lim = Id_ref .* scale;
Iq_ref_lim = Iq_ref .* scale;

log.Iref = [Id_ref; Iq_ref];
log.Iref_mag = sqrt(Id_ref.^2 + Iq_ref.^2);

log.Iref_lim = [Id_ref_lim; Iq_ref_lim];
log.Iref_mag_lim = sqrt(Id_ref_lim.^2 + Iq_ref_lim.^2);

% ============================================================
% Current mag
% ============================================================
log.i1_mag = sqrt(x(1,:).^2 + x(2,:).^2);
log.i2_mag = sqrt(x(3,:).^2 + x(4,:).^2);

% ============================================================
% Converter voltage output
% ============================================================
err_id = Id_ref_lim - i1d;
err_iq = Iq_ref_lim - i1q;

Econv_d = p.Kpc .* err_id + p.Kic .* xi_id ...
            - omega_conv .* p.L1 .* i1q + vPCCd;

Econv_q = p.Kpc .* err_iq + p.Kic .* xi_iq ...
            + omega_conv .* p.L1 .* i1d + vPCCq;


log.Econv = [Econv_d; Econv_q];
log.Econv_mag = sqrt(Econv_d.^2 + Econv_q.^2);

% ============================================================
% Converter power (internal)
% ============================================================
P_conv = Econv_d .* i1d + Econv_q .* i1q;
Q_conv = Econv_q .* i1d - Econv_d .* i1q;

log.P_conv = P_conv;
log.Q_conv = Q_conv;



% ============================================================
% Grid power
% ============================================================
log.P_grid = vg_d .* i2d + vg_q .* i2q;
log.Q_grid = vg_q .* i2d - vg_d .* i2q;

% ============================================================
% Trap filters
% ============================================================
log.P_trap1 = it1d .* vt1d + it1q .* vt1q;
log.Q_trap1 = it1q .* vt1d - it1d .* vt1q;

log.P_trap2 = it2d .* vt2d + it2q .* vt2q;
log.Q_trap2 = it2q .* vt2d - it2d .* vt2q;

% ============================================================
% Angle normalization (plotting)
% ============================================================
shift = delta_g_abs(1);

log.delta_g = delta_g_abs;
log.delta_g_shifted = delta_g_abs - shift;

log.delta_conv = delta_conv;

% ============================================================
% Frequencies
% ============================================================
log.omega_g = omega_g;
log.omega_conv = omega_conv;

end