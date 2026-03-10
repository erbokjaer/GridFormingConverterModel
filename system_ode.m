function dx = system_ode(t,x,p)

% ============================================================
% Unpack states
% ============================================================
i1d   = x(p.state_idx.i1d);
i1q   = x(p.state_idx.i1q);
i2d   = x(p.state_idx.i2d);
i2q   = x(p.state_idx.i2q);
vcd   = x(p.state_idx.vcd);
vcq   = x(p.state_idx.vcq);
delta_g = x(p.state_idx.delta_g);
omega_g = x(p.state_idx.omega_g);
delta_c = x(p.state_idx.delta_c);
omega_c = x(p.state_idx.omega_c);
xi_Q  = x(p.state_idx.xi_Q);
xi_id = x(p.state_idx.xi_id);
xi_iq = x(p.state_idx.xi_iq);
it1d  = x(p.state_idx.it1d);
it1q  = x(p.state_idx.it1q);
vt1d  = x(p.state_idx.vt1d);
vt1q  = x(p.state_idx.vt1q);
it2d  = x(p.state_idx.it2d);
it2q  = x(p.state_idx.it2q);
vt2d  = x(p.state_idx.vt2d);
vt2q  = x(p.state_idx.vt2q);
xi_vd = x(p.state_idx.xi_vd);
xi_vq = x(p.state_idx.xi_vq);


% ============================================================
% Grid voltage events
% ============================================================

vg_mag       = p.vg_mag(t);
vg_phase_rad = p.vg_phase_rad(t);

vg_d = vg_mag * cos(vg_phase_rad);
vg_q = vg_mag * sin(vg_phase_rad);

% ============================================================
% Node voltage
% ============================================================
ilpd = i1d - i2d - it1d - it2d;
ilpq = i1q - i2q - it1q - it2q;

vd = vcd + p.Rlp*ilpd;
vq = vcq + p.Rlp*ilpq;

% ============================================================
% dq transform
% ============================================================
c = cos(delta_c);
s = sin(delta_c);

i1d_c = c*i1d + s*i1q;
i1q_c = -s*i1d + c*i1q;

vd_c = c*vd + s*vq;
vq_c = -s*vd + c*vq;

i2d_c = c*i2d + s*i2q;
i2q_c = -s*i2d + c*i2q;

% ============================================================
% Reactive power control
% ============================================================
V_meas = sqrt(vd_c^2 + vq_c^2);
Q_meas = vq_c*i2d_c - vd_c*i2q_c;

err_v = p.V_ref(t) - V_meas;
err_q = (p.Q_ref(t) - Q_meas) + p.K_vq*err_v;

dxi_Q = err_q;

Vrefd = p.Kp_q*err_q + p.Ki_q*xi_Q;

% ============================================================
% Virtual impedance
% ============================================================

% dVd = Vref - vd_c;
% dVq = -vq_c;
% 
% Id_ref = (dVd*p.R_v_tot + dVq*p.X_v_tot)/p.Z_v2;
% Iq_ref = (dVq*p.R_v_tot - dVd*p.X_v_tot)/p.Z_v2;

e_vd = Vrefd - vd_c;
e_vq = - vq_c;    % same as your previous definition

d_vv_d = p.Kpv*e_vd + p.Kiv*xi_vd;
d_vv_q = p.Kpv*e_vq + p.Kiv*xi_vq;

dxi_vd = e_vd;  % integrator update
dxi_vq = e_vq;

Id_ref = (d_vv_d*p.R_v_tot + d_vv_q*p.X_v_tot)/p.Z_v2;
Iq_ref = (d_vv_q*p.R_v_tot - d_vv_d*p.X_v_tot)/p.Z_v2;

% ============================================================
% Current limit
% ============================================================
I_mag = sqrt(Id_ref^2 + Iq_ref^2);

if I_mag > p.I_max
    scale = p.I_max / I_mag;
    Id_ref = Id_ref * scale;
    Iq_ref = Iq_ref * scale;
end

% ============================================================
% Current controller
% ============================================================
err_id = Id_ref - i1d_c;
err_iq = Iq_ref - i1q_c;

dxi_id = err_id;
dxi_iq = err_iq;

Ecd = p.Kp_c*err_id + p.Ki_c*xi_id ...
      - omega_c*p.L1*i1q_c + vd_c;

Ecq = p.Kp_c*err_iq + p.Ki_c*xi_iq ...
      + omega_c*p.L1*i1d_c + vq_c;

Evd = c*Ecd - s*Ecq;
Evq = s*Ecd + c*Ecq;

% ============================================================
% Electrical dynamics
% ============================================================
di1d=(Evd-vd-p.R1*i1d+omega_g*p.L1*i1q)/p.L1;
di1q=(Evq-vq-p.R1*i1q-omega_g*p.L1*i1d)/p.L1;

di2d=(vd-vg_d-p.R2*i2d+omega_g*p.L2*i2q)/p.L2;
di2q=(vq-vg_q-p.R2*i2q-omega_g*p.L2*i2d)/p.L2;

% ============================================================
% Swing equations
% ============================================================
P_PCC = vd_c*i2d_c + vq_c*i2q_c;

domega_c = (p.P_ref(t) - P_PCC - p.D_conv*(omega_c-p.w_nom))/p.J_conv;
ddelta_c = omega_c - p.w_nom;

P_e = vg_d*i2d + vg_q*i2q;

domega_g=(p.Pm_cont(t) - P_e - p.D_g*(omega_g-p.w_nom))/p.J_g;
ddelta_g=omega_g - p.w_nom;

% ============================================================
% Low-pass filter
% ============================================================
dvcd=(ilpd+omega_g*p.Clp*vcq)/p.Clp;
dvcq=(ilpq-omega_g*p.Clp*vcd)/p.Clp;

% ============================================================
% Trap filters
% ============================================================
dit1d=(vd-vt1d-p.Rt1*it1d+omega_g*p.Lt1*it1q)/p.Lt1;
dit1q=(vq-vt1q-p.Rt1*it1q-omega_g*p.Lt1*it1d)/p.Lt1;

dvt1d=(it1d+omega_g*p.Ct1*vt1q)/p.Ct1;
dvt1q=(it1q-omega_g*p.Ct1*vt1d)/p.Ct1;

dit2d=(vd-vt2d-p.Rt2*it2d+omega_g*p.Lt2*it2q)/p.Lt2;
dit2q=(vq-vt2q-p.Rt2*it2q-omega_g*p.Lt2*it2d)/p.Lt2;

dvt2d=(it2d+omega_g*p.Ct2*vt2q)/p.Ct2;
dvt2q=(it2q-omega_g*p.Ct2*vt2d)/p.Ct2;

% ============================================================
% State vector
% ============================================================
dx = [ ...
    di1d; di1q;            % 1-2: Converter currents
    di2d; di2q;            % 3-4: Grid currents
    dvcd; dvcq;            % 5-6: Capacitor voltage
    ddelta_g; domega_g;    % 7-8: Grid angle & speed
    ddelta_c; domega_c;    % 9-10: Converter angle & speed
    dxi_Q;                 % 11: Reactive power integrator
    dxi_id; dxi_iq;        % 12-13: Current controller integrators
    dit1d; dit1q; dvt1d; dvt1q;  % 14-17: Trap filter 1
    dit2d; dit2q; dvt2d; dvt2q;  % 18-21: Trap filter 2
    dxi_vd; dxi_vq             % 22-23: Virtual impedance integrators
];

end