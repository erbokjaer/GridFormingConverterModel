function log = reconstruct_logs_lin_5_temp(p,t,x)

% ============================================================
% Operating point (clean, no duplicates)
% ============================================================
p.V_PCC0 = p.log0.V_PCC(:,1);
p.Vmag0  = p.log0.V_PCC_mag(1);

p.P0 = p.log0.P_PCC(1);
p.Q0 = p.log0.Q_PCC(1);

p.i10 = [p.x0(p.state_idx.i1d); p.x0(p.state_idx.i1q)];
p.i20 = [p.x0(p.state_idx.i2d); p.x0(p.state_idx.i2q)];

p.vc0  = [p.x0(p.state_idx.vcd); p.x0(p.state_idx.vcq)];

p.Vref0    = p.log0.Vrefd(1);
p.Vref_dq0 = p.log0.Vref(:,1);

p.Econv0 = p.log0.Econv(:,1);

p.Iref0     = p.log0.Iref(:,1);
p.Iref_lim0 = p.log0.Iref_lim(:,1);

p.Pconv0 = p.log0.P_conv;
p.Qconv0 = p.log0.Q_conv;
p.Pgrid0 = p.log0.P_grid;
p.Qgrid0 = p.log0.Q_grid;

p.delta_conv0 = p.x0(p.state_idx.delta_c);
p.omega_conv0 = p.x0(p.state_idx.omega_c);
p.omega_g0    = p.x0(p.state_idx.omega_g);

p.delta_g0 = 0; % removed state → reference frame

Vg0    = p.vg_mag(0);
theta0 = p.vg_phase_rad(0);
p.Vg0  = [Vg0*cos(theta0); Vg0*sin(theta0)];

% ============================================================
% Deviation variables
% ============================================================
dx = x - p.x0;

V_ref  = arrayfun(p.V_ref,t);
P_ref  = arrayfun(p.P_ref,t);
Q_ref  = arrayfun(p.Q_ref,t);
Vg_mag = arrayfun(p.vg_mag,t);
theta_g= arrayfun(p.vg_phase_rad,t);
Pm_cont= arrayfun(p.Pm_cont,t);

du = [V_ref; P_ref; Q_ref] - p.u0;
dd = [Pm_cont; Vg_mag; theta_g] - p.d0;

% ============================================================
% Extract states
% ============================================================
di1d = dx(p.state_idx.i1d,:);
di1q = dx(p.state_idx.i1q,:);
di2d = dx(p.state_idx.i2d,:);
di2q = dx(p.state_idx.i2q,:);
dvcd = dx(p.state_idx.vcd,:);
dvcq = dx(p.state_idx.vcq,:);
ddelta_conv = dx(p.state_idx.delta_c,:);

% ============================================================
% PCC voltage (exact affine)
% ============================================================
dvPCCd = dvcd + p.Rlp*(di1d - di2d ...
    - dx(p.state_idx.it1d,:) - dx(p.state_idx.it2d,:));

dvPCCq = dvcq + p.Rlp*(di1q - di2q ...
    - dx(p.state_idx.it1q,:) - dx(p.state_idx.it2q,:));

log.V_PCC = p.V_PCC0 + [dvPCCd; dvPCCq];

% ============================================================
% Voltage magnitude
% ============================================================
v0d = p.V_PCC0(1);
v0q = p.V_PCC0(2);
V0  = p.Vmag0;

dV = (v0d/V0)*dvPCCd + (v0q/V0)*dvPCCq;
log.V_PCC_mag = V0 + dV;

% ============================================================
% Power (PCC)
% ============================================================
i20d = p.i20(1);
i20q = p.i20(2);

dP = v0d*di2d + i20d*dvPCCd + v0q*di2q + i20q*dvPCCq;
dQ = v0q*di2d + i20d*dvPCCq - v0d*di2q - i20q*dvPCCd;

log.P_PCC = p.P0 + dP;
log.Q_PCC = p.Q0 + dQ;

% ============================================================
% Voltage reference
% ============================================================
dVref_cf = ...
    p.Kpq * (du(3,:) - dQ) + ...
    p.Kpq*p.Kvq * (du(1,:) - dV) + ...
    p.Kiq * dx(p.state_idx.xi_Q,:);

log.Vrefd = p.Vref0 + dVref_cf;

% dq rotation
c0 = cos(p.delta_conv0);
s0 = sin(p.delta_conv0);

dVref_d = c0*dVref_cf - s0*p.Vref0*ddelta_conv;
dVref_q = s0*dVref_cf + c0*p.Vref0*ddelta_conv;

log.Vref = p.Vref_dq0 + [dVref_d; dVref_q];

% ============================================================
% Current reference (linearized, no limiter dynamics)
% ============================================================
dxi_vd = dx(p.state_idx.xi_vd,:);
dxi_vq = dx(p.state_idx.xi_vq,:);

ddv_d = p.Kpv .* dVref_d + p.Kiv .* dxi_vd;
ddv_q = p.Kpv .* dVref_q + p.Kiv .* dxi_vq;


Zv2 = sqrt(p.Rv^2 + p.Xv^2);

dIref_d = (p.Rv .* ddv_d + p.Xv .* ddv_q) ./ Zv2;
dIref_q = (p.Rv .* ddv_q - p.Xv .* ddv_d) ./ Zv2;

log.Iref_lim = p.Iref_lim0 + [dIref_d; dIref_q];

Iref0 = p.Iref0;

Iref_mag0 = norm(Iref0);

dIref_mag = (Iref0(1)*dIref_d + Iref0(2)*dIref_q) / Iref_mag0;

log.Iref_mag_lim = Iref_mag0 + dIref_mag;



% ============================================================
% Current magnitudes
% ============================================================
i10 = p.i10;
i20 = p.i20;

i1_mag0 = norm(i10);
i2_mag0 = norm(i20);

log.i1_mag = i1_mag0 + (i10(1)*di1d + i10(2)*di1q) ./ i1_mag0;
log.i2_mag = i2_mag0 + (i20(1)*di2d + i20(2)*di2q) ./ i2_mag0;

% ============================================================
% Grid voltage (delta_g = 0 frame)
% ============================================================
Vg0 = norm(p.Vg0);
theta0 = 0;

c0 = cos(theta0);
s0 = sin(theta0);

dVg = Vg_mag - Vg0;
dth = theta_g;

dv_gd = c0 .* dVg - Vg0 .* s0 .* dth;
dv_gq = s0 .* dVg + Vg0 .* c0 .* dth;

% ============================================================
% Converter voltage
% ============================================================
dxi_id = dx(p.state_idx.xi_id,:);
dxi_iq = dx(p.state_idx.xi_iq,:);
domega = dx(p.state_idx.omega_c,:);

omega0 = p.omega_conv0;

dEconv_d = ...
    p.Kpc .* (dIref_d - di1d) + ...
    p.Kic .* dxi_id - ...
    p.L1 .* (omega0 .* di1q + i10(2) .* domega) + ...
    dvPCCd;

dEconv_q = ...
    p.Kpc .* (dIref_q - di1q) + ...
    p.Kic .* dxi_iq + ...
    p.L1 .* (omega0 .* di1d + i10(1) .* domega) + ...
    dvPCCq;

log.Econv = p.Econv0 + [dEconv_d; dEconv_q];

E0 = p.Econv0;

E_mag0 = norm(E0);

dE_mag = (E0(1)*dEconv_d + E0(2)*dEconv_q) / E_mag0;

log.Econv_mag = E_mag0 + dE_mag;

% ============================================================
% Converter power
% ============================================================
dP_conv = ...
    p.Econv0(1)*di1d + i10(1)*dEconv_d + ...
    p.Econv0(2)*di1q + i10(2)*dEconv_q;

dQ_conv = ...
    p.Econv0(2)*di1d + i10(1)*dEconv_q - ...
    p.Econv0(1)*di1q - i10(2)*dEconv_d;

log.P_conv = p.Pconv0 + dP_conv;
log.Q_conv = p.Qconv0 + dQ_conv;

% ============================================================
% Grid power
% ============================================================
dP_grid = ...
    p.Vg0(1)*di2d + i20(1)*dv_gd + ...
    p.Vg0(2)*di2q + i20(2)*dv_gq;

dQ_grid = ...
    p.Vg0(2)*di2d + i20(1)*dv_gq - ...
    p.Vg0(1)*di2q - i20(2)*dv_gd;

log.P_grid = p.Pgrid0 + dP_grid;
log.Q_grid = p.Qgrid0 + dQ_grid;

% ============================================================
% Angles & frequencies
% ============================================================
log.delta_g = theta_g; % now purely disturbance
log.delta_conv = p.delta_conv0 + ddelta_conv;

log.omega_g = p.omega_g0 + dx(p.state_idx.omega_g,:);
log.omega_conv = p.omega_conv0 + dx(p.state_idx.omega_c,:);

end