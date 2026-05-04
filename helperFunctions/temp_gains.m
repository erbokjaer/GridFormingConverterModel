
% Map legacy parameter names to the new names used above, so both sets are available

% ============================================================
% Legacy names (from old code) - keep as aliases of the new parameters
% ============================================================
J_g = Jg;        % legacy: generator inertia
D_g = Dg;        % legacy: generator damping

J_conv = Jconv;  % legacy: converter inertia
D_conv = Dconv;  % legacy: converter damping

K_vq = Kvq;      % legacy: reactive voltage gain
Kp_q = Kpq;      % legacy: reactive proportional gain
Ki_q = Kiq;      % legacy: reactive integral gain

% Virtual admittance legacy values: define only if not already set above
% If R_v, X_v, Z_v2 were intended to be set elsewhere, keep them consistent:
if ~exist('R_v','var'); R_v = 0.15; end
if ~exist('X_v','var'); X_v = 0.03; end
Z_v2 = sqrt(R_v^2 + X_v^2);

Kpv = Kpv;       % keep same name (already defined above)
tau_v = tau_v;
Kiv = Kiv;

Kp_c = Kpc;      % legacy: current control proportional
Ki_c = Kic;      % legacy: current control integral

