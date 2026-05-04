function dx = gem_system(t, x, u, d, p)
% GEM nonlinear state-space model
%
% x : state vector (15x1)
% u : [E_conv_d; E_conv_q]
% d : [P_m_cont; Vg_mag; theta_g]
% p : parameter struct

% ============================================================
% Unpack states
% ============================================================
i1d   = x(1);
i1q   = x(2);
i2d   = x(3);
i2q   = x(4);
vcd   = x(5);
vcq   = x(6);

it1d  = x(7);
it1q  = x(8);
vt1d  = x(9);
vt1q  = x(10);

it2d  = x(11);
it2q  = x(12);
vt2d  = x(13);
vt2q  = x(14);

wg    = x(15);

% Inputs
Econv_d = u(1);
Econv_q = u(2);

Pm      = d(1);
Vg_mag  = d(2);
% Vg_mag  = t;
% if t > d(2)
%     Vg_mag  = d(2);
% end

theta_g = d(3);

% ============================================================
% PARAMETERS
% ============================================================
R1   = p.R1;
R2   = p.R2;
Rlp  = p.Rlp;

L1   = p.L1;
L2   = p.L2;

Lt1  = p.Lt1;
Lt2  = p.Lt2;

Ct1  = p.Ct1;
Ct2  = p.Ct2;

Dg   = p.Dg;
Jg   = p.Jg;
% if t > 100
%     1;
% end


% ============================================================
% STATE MATRIX PART (A x)
% ============================================================
A = p.Agem;

Ax = A * x;




% ============================================================
% NONLINEAR TERMS f_nl(x)
% ============================================================
fnl = wg.*p.J*x;

% ============================================================
% g_nl(x) TERMS
% ============================================================
gnl = zeros(16,1);

gnl(3) = -(1/L2) * (Vg_mag) * cos(theta_g);
gnl(4) = -(1/L2) * (Vg_mag) * sin(theta_g);

gnl(15) = (1/Jg) * ( ...
    (Vg_mag)*cos(theta_g)*i2d ...
    + (Vg_mag)*sin(theta_g)*i2q );

% 
% gnl = zeros(15,1);
% 
% gnl(3) = -(1/L2) * ((1 + Vg_mag) * cos(theta_g) - 1) - 1/p.L2;
% gnl(4) = -(1/L2) * (1 + Vg_mag) * sin(theta_g);
% 
% gnl(15) = (1/Jg) * ( ...
%     (1 + Vg_mag)*cos(theta_g)*i2d ...
%     + (1 + Vg_mag)*sin(theta_g)*i2q ) - i2d/Jg;
% 





% ============================================================
% INPUT MATRICES
% ============================================================
B = p.Bgem;
E = p.Egem;

temp = zeros(16,1);
% delta = -x(16);
% u = [cos(delta); sin(delta)];
% 
% temp(16) = 1*wg;

Bu = B * (u);
Ed = E * (d);

% ============================================================
% FINAL DYNAMICS
% ============================================================

dx = Ax + fnl + gnl + Bu + Ed + temp;
% dx(15) = 0;

end