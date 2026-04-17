function log = reconstruct_logs_lin(p,t,x)

log = struct();

% ============================================================
% Deviation variables
% ============================================================
dx = x - p.x0;

% Inputs
V_ref = arrayfun(p.V_ref,t);
Q_ref = arrayfun(p.Q_ref,t);
Vg_mag  = arrayfun(p.vg_mag,t);
theta_g = arrayfun(p.vg_phase_rad,t);
Pm_cont = arrayfun(p.Pm_cont,t);
P_ref = arrayfun(p.P_ref,t);

du = [V_ref; P_ref;Q_ref] - p.u0;
dd = [Pm_cont;Vg_mag; theta_g] - p.d0;

% ============================================================
% Extract deviations (example subset)
% ============================================================
di1d = dx(p.state_idx.i1d,:);
di1q = dx(p.state_idx.i1q,:);
di2d = dx(p.state_idx.i2d,:);
di2q = dx(p.state_idx.i2q,:);
dvcd = dx(p.state_idx.vcd,:);
dvcq = dx(p.state_idx.vcq,:);
ddelta_conv = dx(p.state_idx.delta_c,:);

% ============================================================
% PCC voltage (linearized)
% vPCC = vcd + Rlp*(i1 - i2 - it1 - it2)
% already affine → exact linear form
% ============================================================
dvPCCd = dvcd + p.Rlp*( ...
    di1d ...
  - di2d ...
  - dx(p.state_idx.it1d,:) ...
  - dx(p.state_idx.it2d,:) );

dvPCCq = dvcq + p.Rlp*( ...
    di1q ...
  - di2q ...
  - dx(p.state_idx.it1q,:) ...
  - dx(p.state_idx.it2q,:) );

log.V_PCC = p.V_PCC0 + [dvPCCd; dvPCCq];

% ============================================================
% Voltage magnitude linearization
% sqrt(vd^2 + vq^2)
% ============================================================
v0d = p.V_PCC0(1);
v0q = p.V_PCC0(2);
V0  = p.Vmag0;

dV = (v0d/V0)*dvPCCd + (v0q/V0)*dvPCCq;

log.V_PCC_mag = V0 + dV;

% ============================================================
% Power linearization
% P = vd*i2d + vq*i2q
% ============================================================
i20d = p.i20(1);
i20q = p.i20(2);

dP = v0d*di2d + i20d*dvPCCd + v0q*di2q + i20q*dvPCCq;
dQ = v0q*di2d + i20d*dvPCCq - v0d*di2q - i20q*dvPCCd;

log.P_PCC = p.P0 + dP;
log.Q_PCC = p.Q0 + dQ;

% ============================================================
% Voltage reference (linear)
% ============================================================
dVref_cf = ...
    p.Kpq * (-dQ) + ...
    p.Kpq*p.Kvq * (-dV) + ...
    p.Kiq * dx(p.state_idx.xi_Q,:);

log.Vrefd = p.Vref0 + dVref_cf;

% ============================================================
% dq rotation linearization
% ============================================================
c0 = cos(p.delta_conv0);
s0 = sin(p.delta_conv0);

dVref_d = c0*dVref_cf - s0*p.Vref0*ddelta_conv;
dVref_q = s0*dVref_cf + c0*p.Vref0*ddelta_conv;

log.Vref = p.Vref_dq0 + [dVref_d; dVref_q];

% ============================================================
% Current limit
% ============================================================
% IMPORTANT:
% This is NONLINEAR + NON-SMOOTH → cannot be linearized properly
% Standard approach:
%   Freeze scaling at operating point

% operating point is unsaturated
% Voltage error
dVref = dVref_cf;   % SAME signal used in Vrefd

% PI dynamics contribution (xi_v states must exist)
dxi_vd = dx(p.state_idx.xi_vd,:);
dxi_vq = dx(p.state_idx.xi_vq,:);

% virtual voltage dynamics
ddv_d = p.Kpv .* dVref + p.Kiv .* dxi_vd;
ddv_q = p.Kpv .* dVref + p.Kiv .* dxi_vq;

% impedance mapping (linear transform)
Zv2 = sqrt(p.Rv^2 + p.Xv^2);

dIref_d = (p.Rv .* ddv_d + p.Xv .* ddv_q) ./ Zv2;
dIref_q = (p.Rv .* ddv_q - p.Xv .* ddv_d) ./ Zv2;

% limiter (frozen or removed)
scale0 = 1;

Id_ref_lim = scale0 * dIref_d;
Iq_ref_lim = scale0 * dIref_q;

log.Iref_lim = p.Iref_lim0 + [Id_ref_lim; Iq_ref_lim];


Iref0 = p.Iref0;

Iref_mag = sqrt(Iref0(1,:).^2 + Iref0(2,:).^2);

dIref_mag = (Iref0(1,:).*dIref_d + Iref0(2,:).*dIref_q) ./ Iref_mag;

log.Iref_mag_lim = Iref_mag + dIref_mag;

i10 = p.i10;
i20 = p.i20;

i1_mag0 = sqrt(i10(1,:).^2 + i10(2,:).^2);
i2_mag0 = sqrt(i20(1,:).^2 + i20(2,:).^2);

log.i1_mag = i1_mag0 + ...
    (i10(1,:).*di1d + i10(2,:).*di1q) ./ i1_mag0;

log.i2_mag = i2_mag0 + ...
    (i20(1,:).*di2d + i20(2,:).*di2q) ./ i2_mag0;


% ============================================================
% Converter voltage (linearized)
% ============================================================
% 
% dEconv_d =...
%     p.Kpc*(dIref_d - di1d) + ...
%     p.Kic*dx(p.state_idx.xi_id,:) + ...
%     (-p.L1)*(p.omega_conv0*di1q + p.i10(2,:).*dx(p.state_idx.omega_c,:)) + ...
%     dvPCCd;
% 
% dEconv_q =...
%     p.Kpc*(dIref_q - di1q) + ...
%     p.Kic*dx(p.state_idx.xi_iq,:) + ...
%     ( p.L1)*(p.omega_conv0*di1d + p.i10(1,:).*dx(p.state_idx.omega_c,:)) + ...
%     dvPCCq;
% 
% log.Econv = p.Econv0 + [dEconv_d; dEconv_q];
% 

% 
% Vg0     = p.d0(2);
% theta0  = p.d0(3);

% dVg = Vg_mag - sqrt(p.Vg0(1)^2 - p.Vg0(2)^2);
% dth = theta_g - p.theta_g0;
% dv_gd = c0 .* dVg - Vg0 .* s0 .* dth;
% 
% dv_gq = s0 .* dVg + Vg0 .* c0 .* dth;

dx = x - p.x0;



Vg0 = sqrt(p.Vg0(1)^2 + p.Vg0(2)^2);
p. Vg_mag0 = Vg0;
theta0 = p.theta_g0;


c0 = cos(theta0);
s0 = sin(theta0);


dVg   = Vg_mag - p.Vg_mag0;
dth   = theta_g - p.theta_g0;


dv_gd = c0 .* dVg - Vg0 .* s0 .* dth;
dv_gq = s0 .* dVg + Vg0 .* c0 .* dth;

% vg_d = p.vg0(1,:) + dv_gd;
% vg_q = p.vg0(2,:) + dv_gq;




di1d = dx(p.state_idx.i1d,:);
di1q = dx(p.state_idx.i1q,:);
di2d = dx(p.state_idx.i2d,:);
di2q = dx(p.state_idx.i2q,:);

dit1d = dx(p.state_idx.it1d,:);
dit1q = dx(p.state_idx.it1q,:);
dit2d = dx(p.state_idx.it2d,:);
dit2q = dx(p.state_idx.it2q,:);

% ============================================================
% Converter voltage (linearized correctly)
% ============================================================

% state deviations
di1d = dx(p.state_idx.i1d,:);
di1q = dx(p.state_idx.i1q,:);
dxi_id = dx(p.state_idx.xi_id,:);
dxi_iq = dx(p.state_idx.xi_iq,:);
domega = dx(p.state_idx.omega_c,:);

% you MUST already have these from earlier linearization
% (PCC voltage deviations)
% dvPCCd, dvPCCq

% current reference deviations (must come from your earlier section)
dId_ref = dIref_d;
dIq_ref = dIref_q;

% operating point currents
i10 = p.i10;
omega0 = p.omega_conv0;

% d-axis
dEconv_d = ...
    p.Kpc .* (dId_ref - di1d) + ...
    p.Kic .* dxi_id + ...
    (-p.L1) .* (omega0 .* di1q + i10(2,:) .* domega) + ...
    dvPCCd;

% q-axis
dEconv_q = ...
    p.Kpc .* (dIq_ref - di1q) + ...
    p.Kic .* dxi_iq + ...
    ( p.L1) .* (omega0 .* di1d + i10(1,:) .* domega) + ...
    dvPCCq;

log.Econv = p.Econv0 + [dEconv_d; dEconv_q];

log.Econv_mag = sqrt(p.Econv0(1,:).^2 + p.Econv0(2,:).^2) + ...
    (p.Econv0(1,:).*dEconv_d + p.Econv0(2,:).*dEconv_q) ./ ...
    sqrt(p.Econv0(1,:).^2 + p.Econv0(2,:).^2);


dP_conv = ...
    p.Econv0(1,:) .* di1d + p.i10(1,:) .* dEconv_d + ...
    p.Econv0(2,:) .* di1q + p.i10(2,:) .* dEconv_q;

dQ_conv = ...
    p.Econv0(2,:) .* di1d + p.i10(1,:) .* dEconv_q - ...
    p.Econv0(1,:) .* di1q - p.i10(2,:) .* dEconv_d;

log.P_conv = p.Pconv0 + dP_conv;
log.Q_conv = p.Qconv0 + dQ_conv;




dP_grid = ...
    p.Vg0(1,:) .* di2d + p.i20(1,:) .* dv_gd + ...
    p.Vg0(2,:) .* di2q + p.i20(2,:) .* dv_gq;

dQ_grid = ...
    p.Vg0(2,:) .* di2d + p.i20(1,:) .* dv_gq - ...
    p.Vg0(1,:) .* di2q - p.i20(2,:) .* dv_gd;

log.P_grid = p.Pgrid0 + dP_grid;
log.Q_grid = p.Qgrid0 + dQ_grid;


% dP_trap1 = ...
%     p.it10(1,:) .* dvt1d + p.vt10(1,:) .* dit1d + ...
%     p.it10(2,:) .* dvt1q + p.vt10(2,:) .* dit1q;
% 
% dQ_trap1 = ...
%     p.it10(2,:) .* dvt1d + p.vt10(1,:) .* dit1q - ...
%     p.it10(1,:) .* dvt1q - p.vt10(2,:) .* dit1d;
% 
% log.P_trap1 = p.Ptrap10 + dP_trap1;
% log.Q_trap1 = p.Qtrap10 + dQ_trap1;
% 
% dP_trap2 = ...
%     p.it20(1,:) .* dvt2d + p.vt20(1,:) .* dit2d + ...
%     p.it20(2,:) .* dvt2q + p.vt20(2,:) .* dit2q;
% 
% dQ_trap2 = ...
%     p.it20(2,:) .* dvt2d + p.vt20(1,:) .* dit2q - ...
%     p.it20(1,:) .* dvt2q - p.vt20(2,:) .* dit2d;
% 
% log.P_trap2 = p.Ptrap20 + dP_trap2;
% log.Q_trap2 = p.Qtrap20 + dQ_trap2;

log.delta_g = p.delta_g0 + dx(p.state_idx.delta_g,:);
log.delta_g_shifted = log.delta_g - log.delta_g(:,1);

log.delta_conv = p.delta_conv0 + dx(p.state_idx.delta_c,:);

log.omega_g = p.omega_g0 + dx(p.state_idx.omega_g,:);
log.omega_conv = p.omega_conv0 + dx(p.state_idx.omega_c,:);


end