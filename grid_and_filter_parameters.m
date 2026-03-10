% ============================================================
% Base Values
% ============================================================
S_base = 14e6;
V_base = 690;

w_base = 2*pi*50;
w_nom = 1;

I_base = S_base/(sqrt(3)*V_base);
Z_base = V_base^2/S_base;

L_base = Z_base/w_base;
C_base = 1/(w_base*Z_base);

% ============================================================
% Converter Reactor
% ============================================================
L_reactor_actual = 11.357e-6;
R_reactor_actual = 0.358e-3;

L1 = L_reactor_actual / L_base;
R1 = R_reactor_actual / Z_base;

% ============================================================
% Grid Impedance
% ============================================================
SCR = 2.5;
XR_ratio_grid = 3;

Z2_mag = 1 / SCR;

R2 = Z2_mag / sqrt(1 + XR_ratio_grid^2);
L2 = R2 * XR_ratio_grid;

% ============================================================
% Low Pass Filter
% ============================================================
Clp = 1200e-6 / C_base;
Rlp = 20e-3 / Z_base;

% ============================================================
% Trap Filters
% ============================================================
% % Branch 1: Trap Filter fs (2500 Hz)
Ct1 = 400e-6 / C_base;
Lt1 = 10.3e-6 / L_base;
Rt1 = 1.04e-3 / Z_base;

% Trap Filter 2 (e.g., tuned to ~2500 Hz)
Lt2 = 5.05e-6 / L_base; 
Ct2 = 200e-6 / C_base;
Rt2 = 2.0e-3 / Z_base;