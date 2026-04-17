function dx = system_ode_lin(t, x, p)

% ============================================================
% Inputs
% ============================================================
V_ref = p.V_ref(t);
P_ref = p.P_ref(t); 
Q_ref = p.Q_ref(t);

Pm_cont = p.Pm_cont(t);
Vg_mag  = p.vg_mag(t);
theta_g = p.vg_phase_rad(t); 

u = [V_ref;P_ref;Q_ref] - p.u0;
d = [Pm_cont;Vg_mag;theta_g]- p.d0;

dx = p.Afull * (x - p.x0) + p.Bfull * u + p.Efull * d;%  + p.f_offset;

end