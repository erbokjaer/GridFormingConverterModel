function dx = system_ode_split(t, x, p)

% ============================================================
% Inputs
% ============================================================
V_ref = p.V_ref(t);
P_ref = p.P_ref(t); 
Q_ref = p.Q_ref(t);

Pm_cont = p.Pm_cont(t);
Vg_mag  = p.vg_mag(t);
theta_g = p.vg_phase_rad(t); 

u = [V_ref;P_ref;Q_ref];
d = [Pm_cont;Vg_mag;theta_g];

% ============================================================
% State unpacking
% ============================================================
i1d = x(1);  i1q = x(2);
i2d = x(3);  i2q = x(4);
vcd = x(5);  vcq = x(6);
delta_g = x(7);  omega_g = x(8);

it1d = x(9);  it1q = x(10);
vt1d = x(11); vt1q = x(12);

it2d = x(13); it2q = x(14);
vt2d = x(15); vt2q = x(16);

delta_conv = x(17); omega_conv = x(18);

xi_Q  = x(19);
xi_vd = x(20); xi_vq = x(21);
xi_id = x(22); xi_iq = x(23);

% ============================================================
% PCC voltage
% ============================================================
vPCC_d = vcd + p.Rlp*(i1d - i2d - it1d - it2d);
vPCC_q = vcq + p.Rlp*(i1q - i2q - it1q - it2q);

Vmag = sqrt(vPCC_d^2 + vPCC_q^2);

% ============================================================
% Voltage reference (filtered)
% ============================================================
Vref_cf = ...
    - p.Kpq*(vPCC_q*i2d - vPCC_d*i2q) ...
    - p.Kpq*p.Kvq*Vmag ...
    + p.Kiq*xi_Q;


% ============================================================
% Precompute constants
% ============================================================
Zv2 = sqrt(p.Rv^2 + p.Xv^2);

c = cos(delta_conv);
s = sin(delta_conv);

alpha = (p.Rv*p.Kpv*c + p.Xv*p.Kpv*s)/Zv2;
beta  = (p.Rv*p.Kpv*s - p.Xv*p.Kpv*c)/Zv2;

alpha_c = p.Kpc * alpha;
beta_c  = p.Kpc * beta;

% ============================================================
% Nonlinear f_nl(x)
% ============================================================
fnl = zeros(23,1);

fnl(1) = (alpha_c/p.L1)*Vref_cf - omega_conv*i1q + omega_g*i1q ;
fnl(2) = (beta_c/p.L1)*Vref_cf + omega_conv*i1d - omega_g*i1d ;

fnl(3) = omega_g*i2q;
fnl(4) = -omega_g*i2d;

fnl(5) = omega_g*vcq;
fnl(6) = -omega_g*vcd;

fnl(9)  = omega_g*it1q;
fnl(10) = -omega_g*it1d;

fnl(11) = omega_g*vt1q;
fnl(12) = -omega_g*vt1d;

fnl(13) = omega_g*it2q;
fnl(14) = -omega_g*it2d;

fnl(15) = omega_g*vt2q;
fnl(16) = -omega_g*vt2d;

% % Mechanical dynamics (converter)
fnl(18) = -(1/p.Jconv)*(vPCC_d*i2d + vPCC_q*i2q);

% Q-dynamics
fnl(19) = -(vPCC_q*i2d - vPCC_d*i2q) - p.Kvq*Vmag;

% Voltage controller states
fnl(20) = c * Vref_cf;
fnl(21) = s * Vref_cf;

fnl(22) = alpha * Vref_cf;
fnl(23) = beta  * Vref_cf;

% ============================================================
% g_Q * Q_ref
% ============================================================
gQ = zeros(23,1);

gQ(1) = (alpha_c/p.L1)*p.Kpq;
gQ(2) = (beta_c/p.L1)*p.Kpq;

gQ(20) = c*p.Kpq;
gQ(21) = s*p.Kpq;

gQ(22) = alpha*p.Kpq;
gQ(23) = beta *p.Kpq;

% ============================================================
% g_V * V_ref
% ============================================================
gV = zeros(23,1);

gV(1) = (alpha_c/p.L1)*p.Kpq*p.Kvq;
gV(2) = (beta_c/p.L1)*p.Kpq*p.Kvq;

gV(20) = c*p.Kpq*p.Kvq;
gV(21) = s*p.Kpq*p.Kvq;

gV(22) = alpha*p.Kpq*p.Kvq;
gV(23) = beta *p.Kpq*p.Kvq;

% ============================================================
% Grid disturbance term
% ============================================================
gg = zeros(23,1);

cg = cos(delta_g);
sg = sin(delta_g);

gg(3) = -(1/p.L2)*cg;
gg(4) = -(1/p.L2)*sg;

gg(8) = -(1/p.Jg)*(-cg*i2d - sg*i2q);


% ============================================================
% Total dynamics
% ============================================================
dx = p.A*x + p.B*u + p.E*d + p.f_offset + fnl + gQ*Q_ref + gV*V_ref + gg*Vg_mag;


end