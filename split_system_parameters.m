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
% Rlp = 100e-3 / Z_base;

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


Rv = 0.15;
Xv = 0.03;
Zv2 = sqrt(Rv^2 + Xv^2);

A = zeros(23,23);

% ============================================================
% Row 1
% ============================================================
A(1,1) = (-R1 - Kpc - (Kpc*Rv*Kpv*Rlp)/Zv2)/L1;
A(1,2) = -(Kpc*Xv*Kpv*Rlp)/(L1*Zv2);
A(1,3) =  (Kpc*Rv*Kpv*Rlp)/(L1*Zv2);
A(1,4) =  (Kpc*Xv*Kpv*Rlp)/(L1*Zv2);
A(1,5) = -(Kpc*Rv*Kpv)/(L1*Zv2);
A(1,6) = -(Kpc*Xv*Kpv)/(L1*Zv2);
A(1,9)  = (Kpc*Rv*Kpv*Rlp)/(L1*Zv2);
A(1,10) = (Kpc*Xv*Kpv*Rlp)/(L1*Zv2);
A(1,13) = (Kpc*Rv*Kpv*Rlp)/(L1*Zv2);
A(1,14) = (Kpc*Xv*Kpv*Rlp)/(L1*Zv2);
A(1,20) = (Kpc*Rv*Kiv)/(L1*Zv2);
A(1,21) = (Kpc*Xv*Kiv)/(L1*Zv2);
A(1,22) = Kic/L1;

% ============================================================
% Row 2
% ============================================================
A(2,1) = (Kpc*Xv*Kpv*Rlp)/(L1*Zv2);
A(2,2) = (-R1 - Kpc - (Kpc*Rv*Kpv*Rlp)/Zv2)/L1;
A(2,3) =  -(Kpc*Xv*Kpv*Rlp)/(L1*Zv2);
A(2,4) =  (Kpc*Rv*Kpv*Rlp)/(L1*Zv2);
A(2,5) =  (Kpc*Xv*Kpv)/(L1*Zv2);
A(2,6) = -(Kpc*Rv*Kpv)/(L1*Zv2);
A(2,9)  = -(Kpc*Xv*Kpv*Rlp)/(L1*Zv2);
A(2,10) = (Kpc*Rv*Kpv*Rlp)/(L1*Zv2);
A(2,13) = -(Kpc*Xv*Kpv*Rlp)/(L1*Zv2);
A(2,14) = (Kpc*Rv*Kpv*Rlp)/(L1*Zv2);
A(2,20) = -(Kpc*Xv*Kiv)/(L1*Zv2);
A(2,21) =  (Kpc*Rv*Kiv)/(L1*Zv2);
A(2,23) = Kic/L1;

% ============================================================
% Row 3–6 (grid + capacitor)
% ============================================================
A(3,1) = Rlp/L2;
A(3,3) = -(R2 + Rlp)/L2;
A(3,5) = 1/L2;
A(3,9) = -Rlp/L2;
A(3,13)= -Rlp/L2;

A(4,2) = Rlp/L2;
A(4,4) = -(R2 + Rlp)/L2;
A(4,6) = 1/L2;
A(4,10)= -Rlp/L2;
A(4,14)= -Rlp/L2;

A(5,1) = 1/Clp;
A(5,3) = -1/Clp;
A(5,9) = -1/Clp;
A(5,13)= -1/Clp;

A(6,2) = 1/Clp;
A(6,4) = -1/Clp;
A(6,10)= -1/Clp;
A(6,14)= -1/Clp;

% ============================================================
% Grid angle dynamics
% ============================================================
A(7,8) = 1;
A(8,8) = -Dg/Jg;

% ============================================================
% Trap filter 1
% ============================================================
A(9,1) = Rlp/Lt1;
A(9,3) = -Rlp/Lt1;
A(9,5) = 1/Lt1;
A(9,9) = -(Rlp + Rt1)/Lt1;
A(9,11)= -1/Lt1;
A(9,13)= -Rlp/Lt1;

A(10,2) = Rlp/Lt1;
A(10,4) = -Rlp/Lt1;
A(10,6) = 1/Lt1;
A(10,10)= -(Rlp + Rt1)/Lt1;
A(10,12)= -1/Lt1;
A(10,14)= -Rlp/Lt1;

A(11,9) = 1/Ct1;
A(12,10)= 1/Ct1;

% ============================================================
% Trap filter 2
% ============================================================
A(13,1) = Rlp/Lt2;
A(13,3) = -Rlp/Lt2;
A(13,5) = 1/Lt2;
A(13,9) = -Rlp/Lt2;
A(13,13)= -(Rlp + Rt2)/Lt2;
A(13,15)= -1/Lt2;

A(14,2) = Rlp/Lt2;
A(14,4) = -Rlp/Lt2;
A(14,6) = 1/Lt2;
A(14,10)= -Rlp/Lt2;
A(14,14)= -(Rlp + Rt2)/Lt2;
A(14,16)= -1/Lt2;

A(15,13)= 1/Ct2;
A(16,14)= 1/Ct2;

% ============================================================
% Converter angle
% ============================================================
A(17,18) = 1;
A(18,18) = -Dconv/Jconv;

% ============================================================
% Algebraic-like rows (controllers)
% ============================================================
A(20,1) = -Rlp;
A(20,3) = Rlp;
A(20,5) = -1;
A(20,9) = Rlp;
A(20,13)= Rlp;

A(21,2) = -Rlp;
A(21,4) = Rlp;
A(21,6) = -1;
A(21,10)= Rlp;
A(21,14)= Rlp;

% ============================================================
% Voltage controller internal states
% ============================================================
A(22,1) = -1-(Rv*Kpv*Rlp)/Zv2;
A(22,2) = -(Xv*Kpv*Rlp)/Zv2;
A(22,3) =  (Rv*Kpv*Rlp)/Zv2;
A(22,4) =  (Xv*Kpv*Rlp)/Zv2;
A(22,5) = -(Rv*Kpv)/Zv2;
A(22,6) = -(Xv*Kpv)/Zv2;
A(22,9) =  (Rv*Kpv*Rlp)/Zv2;
A(22,10)=  (Xv*Kpv*Rlp)/Zv2;
A(22,13)=  (Rv*Kpv*Rlp)/Zv2;
A(22,14)=  (Xv*Kpv*Rlp)/Zv2;
A(22,20)=  (Rv*Kiv)/Zv2;
A(22,21)=  (Xv*Kiv)/Zv2;

A(23,1) = (Xv*Kpv*Rlp)/Zv2;
A(23,2) = -1-(Rv*Kpv*Rlp)/Zv2;
A(23,3) = -(Xv*Kpv*Rlp)/Zv2;
A(23,4) =  (Rv*Kpv*Rlp)/Zv2;
A(23,5) = (Xv*Kpv)/Zv2;
A(23,6) = -(Rv*Kpv)/Zv2;
A(23,9) =  -(Xv*Kpv*Rlp)/Zv2;
A(23,10)=  (Rv*Kpv*Rlp)/Zv2;
A(23,13)=  -(Xv*Kpv*Rlp)/Zv2;
A(23,14)=  (Rv*Kpv*Rlp)/Zv2;
A(23,20)= -(Xv*Kiv)/Zv2;
A(23,21)=  (Rv*Kiv)/Zv2;


B = zeros(23,3);

% Inputs: [V_ref; P_ref; Q_ref]

B(18,2) = 1/Jconv;  
B(19,3) = 1;           
B(19,1) = Kvq;       


E = zeros(23,3);

% Disturbances: [Pm_cont; Vg_mag; theta_g]

E(8,1) = -1/Jg;   % mechanical input to grid frequency


f_offset = zeros(23,1);

% ============================================================
% Grid angle dynamics offset
% ============================================================
f_offset(7) = -w_nom;

% Swing equation offset (grid)
f_offset(8) = (Dg * w_nom) / Jg;

% ============================================================
% Converter angle dynamics offset
% ============================================================
f_offset(17) = -w_nom;

% Swing equation offset (converter)
f_offset(18) = (Dconv * w_nom) / Jconv;