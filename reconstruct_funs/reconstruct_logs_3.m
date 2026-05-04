function log = reconstruct_logs_3(p,t,x)

log = struct();

% ============================================================
% State unpacking
% ============================================================
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

V_ref = arrayfun(p.V_ref,t);
Q_ref = arrayfun(p.Q_ref,t);
vg_mag       = arrayfun(p.vg_mag,t);
vg_phase_rad = arrayfun(p.vg_phase_rad,t);
delta_g = x(p.state_idx.delta_g,:) + vg_phase_rad;

% ============================================================
% Grid voltage events
% ============================================================


% vg_d = vg_mag .* cos(vg_phase_rad);
% vg_q = vg_mag .* sin(vg_phase_rad);
vg_d = vg_mag .* cos(delta_g);
vg_q = vg_mag .* sin(delta_g);

log.vg_d = vg_d;
log.vg_q = vg_q;

% ============================================================
% Node voltage (PCC)
% ============================================================
ilpd = i1d - i2d - it1d - it2d;
ilpq = i1q - i2q - it1q - it2q;

vPCCd = vcd + p.Rlp .* ilpd;
vPCCq = vcq + p.Rlp .* ilpq;

log.V_PCC = [vPCCd; vPCCq];
log.V_PCC_mag = sqrt(vPCCd.^2 + vPCCq.^2);


% ============================================================
% Reactive power controller signals
% ============================================================
V_meas = sqrt(vPCCd.^2 + vPCCq.^2);

Q_meas = vPCCq .* i2d - vPCCd .* i2q;
err_v = V_ref - V_meas;
err_q = (Q_ref - Q_meas) + p.K_vq .* err_v;

log.Vrefd = p.Kp_q .* err_q + p.Ki_q .* xi_Q;
 
Vrefd_c = p.Kp_q .* err_q + p.Ki_q .* xi_Q;
Vrefq_c = 0;

c = cos(delta_conv);
s = sin(delta_conv);
Vrefd  =  c.*Vrefd_c - s.*Vrefq_c;
Vrefq = s.*Vrefd_c + c.*Vrefq_c;

% ============================================================
% Virtual impedance current reference
% ============================================================
e_vd = Vrefd - vPCCd;
e_vq = Vrefq - vPCCq;  

d_vv_d = p.Kpv.*e_vd + p.Kiv.*xi_vd;
d_vv_q = p.Kpv.*e_vq + p.Kiv.*xi_vq;


Id_ref = (d_vv_d .* p.R_v + d_vv_q .* p.X_v)./p.Z_v2;
Iq_ref = (d_vv_q .* p.R_v - d_vv_d .* p.X_v)./p.Z_v2;

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

Econv_d = p.Kp_c .* err_id + p.Ki_c .* xi_id ...
            - omega_conv .* p.L1 .* i1q + vPCCd;

Econv_q = p.Kp_c .* err_iq + p.Ki_c .* xi_iq ...
            + omega_conv .* p.L1 .* i1d + vPCCq;


log.Econv = [Econv_d; Econv_q];
log.Econv_mag = sqrt(Econv_d.^2 + Econv_q.^2);

% ============================================================
% Power calculations
% ============================================================
log.P_conv = Econv_d .* i1d + Econv_q .* i1q;
log.Q_conv = Econv_q .* i1d - Econv_d .* i1q;

log.P_PCC = vPCCd .* i2d + vPCCq .* i2q;
log.Q_PCC = vPCCq .* i2d - vPCCd .* i2q;

log.P_grid = vg_d .* i2d + vg_q .* i2q;
log.Q_grid = vg_q .* i2d - vg_d .* i2q;

% Trap filters
log.P_trap1 = it1d .* vt1d + it1q .* vt1q;
log.Q_trap1 = it1q .* vt1d - it1d .* vt1q;

log.P_trap2 = it2d .* vt2d + it2q .* vt2q;
log.Q_trap2 = it2q .* vt2d - it2d .* vt2q;

% ============================================================
% Angle shift for plotting
% ============================================================
log.shift_val = delta_g(1);
log.delta_g_shifted = delta_g - log.shift_val;

% ============================================================
% Synchronization states
% ============================================================
log.delta_g = delta_g;
log.delta_conv = delta_conv;
log.omega_g = omega_g;
log.omega_conv = omega_conv;


end