
% ============================================================
% Symbolic variables
% ============================================================
x  = sym('x',  [23 1], 'real');   % states
u  = sym('u',  [3 1],  'real');   % inputs: [V_ref; P_ref; Q_ref]
d  = sym('d',  [3 1],  'real');   % disturbances: [Pm; Vg; theta_g]

% Unpack states (for readability)
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

% Inputs
V_ref = u(1);
P_ref = u(2);
Q_ref = u(3);

% Disturbances
Pm_cont = d(1);
Vg_mag  = d(2);
theta_g = d(3);

% ============================================================
% Parameters (symbolic)
% ============================================================
p = struct();

params = {'Rlp','Kpq','Kvq','Kiq','Rv','Xv','Kpv','Kpc',...
          'L1','L2','Jconv','Jg'};

for k = 1:length(params)
    p.(params{k}) = sym(params{k}, 'real');
end

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
% Nonlinear dynamics f(x,u,d)
% ============================================================
f = sym(zeros(23,1));

f(1) = (alpha_c/p.L1)*Vref_cf - omega_conv*i1q + omega_g*i1q;
f(2) = (beta_c/p.L1)*Vref_cf + omega_conv*i1d - omega_g*i1d;

f(3) = omega_g*i2q;
f(4) = -omega_g*i2d;

f(5) = omega_g*vcq;
f(6) = -omega_g*vcd;

f(7) = omega_g;                % delta_g dot
f(8) = -(1/p.Jconv)*(vPCC_d*i2d + vPCC_q*i2q);

f(9)  = omega_g*it1q;
f(10) = -omega_g*it1d;

f(11) = omega_g*vt1q;
f(12) = -omega_g*vt1d;

f(13) = omega_g*it2q;
f(14) = -omega_g*it2d;

f(15) = omega_g*vt2q;
f(16) = -omega_g*vt2d;

f(17) = omega_conv;            % delta_conv dot
f(18) = -(1/p.Jconv)*(vPCC_d*i2d + vPCC_q*i2q);

f(19) = -(vPCC_q*i2d - vPCC_d*i2q) - p.Kvq*Vmag;

f(20) = c * Vref_cf;
f(21) = s * Vref_cf;

f(22) = alpha * Vref_cf;
f(23) = beta  * Vref_cf;

% ============================================================
% Add input and disturbance structure
% ============================================================
% gQ * Q_ref
gQ = sym(zeros(23,1));
gQ([1 2]) = [(alpha_c/p.L1)*p.Kpq;
             (beta_c/p.L1)*p.Kpq];

gQ([20 21]) = [c*p.Kpq;
               s*p.Kpq];

gQ([22 23]) = [alpha*p.Kpq;
               beta*p.Kpq];

% gV * V_ref
gV = sym(zeros(23,1));
gV([1 2]) = [(alpha_c/p.L1)*p.Kpq*p.Kvq;
             (beta_c/p.L1)*p.Kpq*p.Kvq];

gV([20 21]) = [c*p.Kpq*p.Kvq;
               s*p.Kpq*p.Kvq];

gV([22 23]) = [alpha*p.Kpq*p.Kvq;
               beta*p.Kpq*p.Kvq];

% disturbance gg * Vg
gg = sym(zeros(23,1));

cg = cos(delta_g);
sg = sin(delta_g);

gg(3) = -(1/p.L2)*cg;
gg(4) = -(1/p.L2)*sg;

gg(8) = -(1/p.Jg)*(-cg*i2d - sg*i2q);

% ============================================================
% Full dynamics
% ============================================================
f_full = f + gQ*Q_ref + gV*V_ref + gg*Vg_mag;

% ============================================================
% Linearization (Jacobian)
% ============================================================
disp('Computing Jacobians...')

Anl = jacobian(f_full, x);
Bnl = jacobian(f_full, u);
Enl = jacobian(f_full, d);

% % ============================================================
% % Convert to MATLAB functions (optional but recommended)
% % ============================================================
matlabFunction(Anl, 'File', 'A_fun', 'Vars', {x,u,d,struct2array(p)}, 'Optimize',true);
matlabFunction(Bnl, 'File', 'B_fun', 'Vars', {x,u,d,struct2array(p)}, 'Optimize',true);
matlabFunction(Enl, 'File', 'E_fun', 'Vars', {x,u,d,struct2array(p)}, 'Optimize',true);