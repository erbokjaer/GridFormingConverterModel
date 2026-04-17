% ============================================================
% Setup
% ============================================================
run setup.m



% run liniarize.m

run grid_and_filter_parameters.m
run temp_gains.m
p = struct;

% grid/filter parameters automatically available from scripts
vars = who;
for k = 1:length(vars)
    p.(vars{k}) = eval(vars{k});
end

use_last_end_x_as_init = 1;
save_last_end_x_as_init = 0;

init_zero = 0;
find_steady_state = 1;
sys_num = 1;
if sys_num == 1
    sys_ode = @system_ode_split;
    reconstruct_fun = @reconstruct_logs_split;
elseif sys_num == 2
    sys_ode = @system_ode_3;
    reconstruct_fun = @reconstruct_logs_3;
elseif sys_num == 3
    sys_ode = @system_ode;
    reconstruct_fun = @reconstruct_logs;
end

% update_plots = [
%     1; % 1
%     0; % 2
%     1; % 3
%     0; % 4
%     0; % 5
%     0; % 6
%     0; % 7
% ];

update_plots = 1*ones(1,7);

% ============================================================
% Simulation horizon
% ============================================================
t_end = 500;
tspan = [0 t_end];

% ============================================================
% References
% ============================================================
% Grid voltage 
Vg.y_prefault  = 1.0;    % normal voltage
Vg.y_fault     = 1;    % voltage during fault
Vg.y_postfault = 1.0;    % voltage after clearing

Vg.t_fault     = 100;      % fault start
Vg.t_apply     = 0.;   % ramp duration into fault

Vg.t_clear     = 101;   % clearing time
Vg.t_recover   = 0.;   % ramp duration back

% Grid phase 
Phg.t_start = 150;
Phg.t_dur   = 0;
Phg.y0      = 0;     % degrees
Phg.y1      = 0;     % degrees

% Voltage reference
V.t_start = 40;
V.t_dur   = 10;
V.y0      = 1.00;
V.y1      = 1.00;

% Reactive power reference
Q.t_start = 20;
Q.t_dur   = 0;
Q.y0      = 0.0;
Q.y1      = 0.0;

% Active power reference
P.t_start = 10;
P.t_dur   = 0;
P.y0      = 0.5;
P.y1      = 0.7;

% Mechanical power
Pm.t_start = 10;
Pm.t_dur   = 0;
Pm.y0      = 0.5;
Pm.y1      = 0.7;

vg_mag = @(t) fault_profile(t,Vg);
vg_phase_rad = @(t) ramp_signal(t, Phg.t_start, Phg.t_dur, deg2rad(Phg.y0), deg2rad(Phg.y1));
V_ref = @(t) ramp_signal(t, V.t_start,  V.t_dur,  V.y0,  V.y1);
Q_ref = @(t) ramp_signal(t, Q.t_start,  Q.t_dur,  Q.y0,  Q.y1);
P_ref     = @(t) ramp_signal(t, P.t_start,  P.t_dur,  P.y0,  P.y1);
Pm_cont   = @(t) ramp_signal(t, Pm.t_start, Pm.t_dur, Pm.y0, Pm.y1);

if init_zero
    vg_mag = @(t) 1;
    vg_phase_rad = @(t) 0;
    V_ref = @(t) 1;
    Q_ref = @(t) 0;
    P_ref = @(t) 0;
    Pm_cont = @(t) 0;
end

% ============================================================
% Faults and limits
% ============================================================
I_max = 10.2;

% ============================================================
% Swing Equation Parameters Grid
% ============================================================
Jg = 1;
Dg = 70;

% ============================================================
% Swing Equation Parameters Converter
% ============================================================
Jconv = 10;
Dconv = 70;

% ============================================================
% Reactive Power Control
% ============================================================
Kvq =  10;      
Kpq = 10; 
tau_q_i = 0.5;
Kiq = Kvq/tau_q_i; 

% ============================================================
% Virtual Admittance
% ============================================================
% Set in parameter file:
% R_v = 0.15;
% X_v = 0.03;
% Z_v2 = sqrt(R_v^2 + X_v^2);

% Virtual impedance PI gains
Kpv = 1;
tau_v = 0.25;
Kiv = Kpv/tau_v;

% ============================================================
% Current Control
% ============================================================
Kpc = 1; 
tau_c_i = 0.1;
Kic = Kpc/tau_c_i; 

% ============================================================
% Load Grid and Filter Parameters
% ============================================================
run split_system_parameters.m


% ============================================================
% Initial state
% ============================================================
% State indices
state_idx = struct( ...
    'i1d', 1, 'i1q', 2, ...         % Converter currents
    'i2d', 3, 'i2q', 4, ...         % Grid currents
    'vcd', 5, 'vcq', 6, ...         % Capacitor voltage
    'delta_g', 7, 'omega_g', 8, ... % Grid generator angle & speed
    'it1d', 9, 'it1q', 10, ...      % Trap filter 1
    'vt1d', 11, 'vt1q', 12, ...
    'it2d', 13, 'it2q', 14, ...     % Trap filter 2
    'vt2d', 15, 'vt2q', 16, ...
    'delta_c', 17,'omega_c', 18,... % Converter angle & speed
    'xi_Q', 19, ...                 % Reactive power integrator
    'xi_vd', 20, 'xi_vq', 21, ...   % Virtual impedance integrators
    'xi_id', 22, 'xi_iq', 23 ...    % Current controller integrators
);

n_states = numel(fieldnames(state_idx));

x0 = zeros(n_states,1);

x0(state_idx.vcd) = 1;
x0(state_idx.omega_g) = w_nom;
x0(state_idx.omega_c) = w_nom;
x0(state_idx.vt1d) = 1;
x0(state_idx.vt2d) = 1;


% ============================================================
% Collect parameters in struct
% ============================================================
p = struct;

% grid/filter parameters automatically available from scripts
vars = who;
for k = 1:length(vars)
    p.(vars{k}) = eval(vars{k});
end


% ============================================================
% Find steady-state initial condition (optional)
% ============================================================
res = inf;
checkpoint_file = 'init_data/system_ode_states.mat';
if use_last_end_x_as_init && exist(checkpoint_file,'file')
    load(checkpoint_file);
    x0 = last_x;
    res = norm(sys_ode(0, x0, p));
end



if find_steady_state && (res > 1e-7)
    f_ss = @(x) sys_ode(0, x, p);
    opts = optimoptions('lsqnonlin', ...
        'Display','off', ...           
        'FunctionTolerance',1e-12, ...
        'StepTolerance',1e-12, ...
        'OptimalityTolerance',1e-12, ...
        'MaxFunctionEvaluations',1e6);
    try
        x_ss = lsqnonlin(@(x) f_ss(x), x0, [], [],opts);
        res = norm(f_ss(x_ss));
        if res < 1e-7
            x0 = x_ss;
        else
            x0 = x_ss;
            warning('Steady-state solver did not converge to sufficient tolerance (residual = %g).', res);
        end
    catch ME
        warning(['Steady-state solver failed: ', ME.message, '. Using original x0.']);
    end
end

x0(state_idx.delta_c) = mod(x0(state_idx.delta_c), 2*pi);


u0 = [V_ref(0); P_ref(0); Q_ref(0)];
d0 = [Pm_cont(0); vg_mag(0); vg_phase_rad(0)];

p.x0 = x0;
p.u0 = u0;
p.d0 = d0;

p_num = [ ...
    Rlp;
    Kpq;
    Kvq;
    Kiq;
    Rv;
    Xv;
    Kpv;
    Kpc;
    L1;
    L2;
    Jconv;
    Jg
]';


dx = sys_ode(0, x0, p);
A_num = A_fun(x0, u0, d0, p_num);
max(real(eig(A_num)));


Afull = A + A_fun(x0, u0, d0, p_num);
Bfull = B + double(B_fun(x0,u0,d0,p_num));
Efull = E + double(E_fun(x0,u0,d0,p_num));
lambda_Afull = eig(Afull);


fig98 = figure(98);
% theme(fig98, 'light')
plot(real(lambda_Afull), imag(lambda_Afull), 'x', LineWidth=2);
grid on;

xlabel('Real Part');
ylabel('Imaginary Part');
title('Eigenvalues of A_{full}');



p.Afull = Afull;
p.Bfull = Bfull;
p.Efull = Efull;

p_vec = struct();

params = {'Rlp','Kpq','Kvq','Kiq','Rv','Xv','Kpv','Kpc',...
          'L1','L2','Jconv','Jg'};

for k = 1:length(params)
    p_vec.(params{k}) = sym(params{k}, 'real');
end



% A_eval = double(subs(Anl, ...
%     [x; u; d; struct2array(p_vec)'], ...
%     [x0; u0; d0; p_num']));
% 
% max(abs(A_eval - A_num));



% ============================================================
% Solver
% ============================================================
% opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-2);
opts = odeset('RelTol',1e-4, 'AbsTol',1e-8, 'MaxStep',1e-1);

% [t,x] = ode23t(@(t,x) sys_ode(t,x,p), tspan, x0, opts);
% [t,x] = ode15s(@(t,x) sys_ode(t,x,p), tspan, x0, opts);
% [t,x] = ode45(@(t,x) sys_ode(t,x,p), tspan, x0, opts);

t_fixed = linspace(tspan(1), tspan(end), 100000);
[t1,x1] = ode23t(@(t,x) system_ode_split(t,x,p), t_fixed, x0, opts);
% [t2,x2] = ode23t(@(t,x) system_ode(t,x,p), t_fixed, x0, opts);

% [t1,x1] = ode23t(@(t,x) system_ode_3(t,x,p), t_fixed, x0, opts);

% [t1,x1] = ode23t(@(t,x) system_ode(t,x,p), t_fixed, x0, opts);
[t2,x2] = ode23t(@(t,x) system_ode_lin(t,x,p), t_fixed, x0, opts);
% [t2,x2] = ode23t(@(t,x) system_ode_4(t,x,p), t_fixed, x0, opts);
t = t1.';         
x1 = x1.';        
x2 = x2.';

log1 = reconstruct_logs_split(p,t,x1);
% STATE_NUM = 1:2;
% x2(STATE_NUM,:) = x1(STATE_NUM,:);

log0 = reconstruct_logs_split(p, 0, p.x0);
% PCC voltage
p.V_PCC0 = log0.V_PCC(:,1);
p.Vmag0  = log0.V_PCC_mag(1);

% Power
p.P0 = log0.P_PCC(1);
p.Q0 = log0.Q_PCC(1);

% Current
p.i10 = [p.x0(p.state_idx.i1d); p.x0(p.state_idx.i1q)];
p.i20 = [p.x0(p.state_idx.i2d); p.x0(p.state_idx.i2q)];

% Voltage reference
p.Vref0     = log0.Vrefd(1);
p.Vref_dq0  = log0.Vref(:,1);

% Converter voltage
p.Econv0 = log0.Econv(:,1);

% Current reference (before and after limiter)
p.Iref0      = log0.Iref(:,1);
p.Iref_lim0  = log0.Iref_lim(:,1);

% Limiter scaling
p.scale0 = norm(p.Iref_lim0) / max(norm(p.Iref0),1e-12);

% Angles / frequencies
p.delta_conv0 = p.x0(p.state_idx.delta_c);
p.omega_conv0 = p.x0(p.state_idx.omega_c);

p.Econv0 = log0.Econv;
% =========================
% STATES
% =========================
p.i10  = [x0(p.state_idx.i1d); x0(p.state_idx.i1q)];
p.i20  = [x0(p.state_idx.i2d); x0(p.state_idx.i2q)];

p.it10 = [x0(p.state_idx.it1d); x0(p.state_idx.it1q)];
p.it20 = [x0(p.state_idx.it2d); x0(p.state_idx.it2q)];

p.vc0  = [x0(p.state_idx.vcd); x0(p.state_idx.vcq)];

p.Pconv0 = log0.P_conv;
p.Qconv0 = log0.Q_conv;
p.Pgrid0 = log0.P_grid;
p.Qgrid0 = log0.Q_grid;

Vg0    = p.vg_mag(0);
p.theta_g0 = p.vg_phase_rad(0);

p.Vg0 = [Vg0*cos(p.theta_g0); Vg0*sin(p.theta_g0)];

% =========================
% FREQUENCIES / ANGLES
% =========================
p.delta_conv0 = x0(p.state_idx.delta_c);
p.omega_conv0 = x0(p.state_idx.omega_c);

p.delta_g0 = x0(p.state_idx.delta_g);
p.omega_g0 = x0(p.state_idx.omega_g);

log2 = reconstruct_logs_lin(p,t,x2);





if save_last_end_x_as_init && not(any(isnan(x1(:, end))))
    last_x = x1(:, end);
    save(checkpoint_file, 'last_x');

end



% ============================================================
% Fig 1: Currents and PCC Voltage
% ============================================================
if update_plots(1)
    figure(1); clf; set(gcf,'WindowStyle','docked');

    subplot(3,1,1)
    plot(t,x1(state_idx.i1d,:),'-', ...
         t,x1(state_idx.i1q,:),'-', ...
         t,x2(state_idx.i1d,:),'--', ...
         t,x2(state_idx.i1q,:),'--','LineWidth',1)
    legend('$i_{1d}$','$i_{1q}$','$i_{1d}^{(2)}$','$i_{1q}^{(2)}$')
    grid on; title('Converter Currents'); axis padded

    subplot(3,1,2)
    plot(t,x1(state_idx.i2d,:),'-', ...
         t,x1(state_idx.i2q,:),'-', ...
         t,x2(state_idx.i2d,:),'--', ...
         t,x2(state_idx.i2q,:),'--','LineWidth',1)
    legend('$i_{2d}$','$i_{2q}$','$i_{2d}^{(2)}$','$i_{2q}^{(2)}$')
    grid on; title('Grid Side Currents'); axis padded

    subplot(3,1,3)
    plot(t,log1.V_PCC(1,:),'-', ...
         t,log1.V_PCC(2,:),'-', ...
         t,log2.V_PCC(1,:),'--', ...
         t,log2.V_PCC(2,:),'--','LineWidth',1)
    legend('$v_d$','$v_q$','$v_d^{(2)}$','$v_q^{(2)}$')
    grid on; title('PCC Voltage'); axis padded
end

% ============================================================
% Fig 2: Active and Reactive Power
% ============================================================

if update_plots(2)
    figure(2); clf; set(gcf,'WindowStyle','docked');

    subplot(2,1,1)
    plot(t,log1.P_conv,'-', t,log2.P_conv,'--', ...
         t,log1.P_grid,'-', t,log2.P_grid,'--', ...
         t,log1.P_PCC,'-', t,log2.P_PCC,'--','LineWidth',1)
    legend('$P_{conv}$','$P_{conv}^{(2)}$', ...
           '$P_{grid}$','$P_{grid}^{(2)}$', ...
           '$P_{PCC}$','$P_{PCC}^{(2)}$')
    grid on; title('Active Power'); axis padded

    subplot(2,1,2)
    plot(t,log1.Q_conv,'-', t,log2.Q_conv,'--', ...
         t,log1.Q_grid,'-', t,log2.Q_grid,'--', ...
         t,log1.Q_PCC,'-', t,log2.Q_PCC,'--','LineWidth',1)
    legend('$Q_{conv}$','$Q_{conv}^{(2)}$', ...
           '$Q_{grid}$','$Q_{grid}^{(2)}$', ...
           '$Q_{PCC}$','$Q_{PCC}^{(2)}$')
    grid on; title('Reactive Power'); axis padded
end

% ============================================================
% Fig 3: Synchronization States
% ============================================================
if update_plots(3)
    figure(3); clf; set(gcf,'WindowStyle','docked');

    subplot(2,1,1)
    plot(t,log1.delta_g_shifted,'-', ...
         t,log2.delta_g_shifted,'--', ...
         t,log1.delta_conv,'-', ...
         t,log2.delta_conv,'--','LineWidth',1.2)
    legend('$\delta_g$','$\delta_g^{(2)}$', ...
           '$\delta_{conv}$','$\delta_{conv}^{(2)}$')
    grid on; axis padded

    subplot(2,1,2)
    plot(t,log1.omega_g,'-', ...
         t,log2.omega_g,'--', ...
         t,log1.omega_conv,'-', ...
         t,log2.omega_conv,'--','LineWidth',1.2)
    legend('$\omega_g$','$\omega_g^{(2)}$', ...
           '$\omega_{conv}$','$\omega_{conv}^{(2)}$')
    grid on; axis padded
    title('Synchronization States'); xlabel('Time (s)')
end

% ============================================================
% Fig 4: Vref and PCC Voltage
% ============================================================
if update_plots(4)
    figure(4); clf; set(gcf,'WindowStyle','docked');

    plot(t,log1.Vrefd,'-', ...
         t,log2.Vrefd,'--', ...
         t,log1.V_PCC_mag,'-', ...
         t,log2.V_PCC_mag,'--','LineWidth',1.2)
    grid on
    title('Voltage Reference and PCC Voltage')
    xlabel('Time (s)'); ylabel('Voltage (pu)')
    legend('$V_{ref,d}$','$V_{ref,d}^{(2)}$', ...
           '$V_{PCC}$','$V_{PCC}^{(2)}$')
    axis padded
end

% ============================================================
% Fig 5: Reference Current Magnitude
% ============================================================
if update_plots(5)
    figure(5); clf; set(gcf,'WindowStyle','docked');

    plot(t,log1.Iref_mag_lim,'-', ...
         t,log2.Iref_mag_lim,'--','LineWidth',1.2)
    grid on
    title('Reference Current Magnitude')
    xlabel('Time (s)'); ylabel('$I_{ref}$ (pu)')
    legend('$I_{ref}$','$I_{ref}^{(2)}$')
    axis padded
end

% ============================================================
% Fig 6: Current Magnitudes
% ============================================================
if update_plots(6)
    figure(6); clf; set(gcf,'WindowStyle','docked');

    plot(t,log1.i1_mag,'-', t,log2.i1_mag,'--', ...
         t,log1.i2_mag,'-', t,log2.i2_mag,'--','LineWidth',1.2)
    grid on
    title('Current Magnitudes')
    xlabel('Time (s)'); ylabel('Current (pu)')
    legend('$|i_1|$','$|i_1|^{(2)}$', ...
           '$|i_2|$','$|i_2|^{(2)}$')
    axis padded
end

% ============================================================
% Fig 7: Converter Internal Voltage
% ============================================================
if update_plots(7)
    figure(7); clf; set(gcf,'WindowStyle','docked');

    subplot(2,1,1)
    plot(t,log1.Econv(1,:),'-', ...
         t,log1.Econv(2,:),'-', ...
         t,log2.Econv(1,:),'--', ...
         t,log2.Econv(2,:),'--','LineWidth',1.2)
    legend('$E_{conv,d}$','$E_{conv,q}$', ...
           '$E_{conv,d}^{(2)}$','$E_{conv,q}^{(2)}$')
    grid on; title('Converter Internal Voltage (d/q)')
    xlabel('Time (s)'); axis padded

    subplot(2,1,2)
    plot(t,log1.Econv_mag,'-', ...
         t,log2.Econv_mag,'--','LineWidth',1.2)
    grid on
    title('Converter Internal Voltage Magnitude')
    xlabel('Time (s)'); ylabel('Voltage (pu)')
    legend('$|E_{conv}|$','$|E_{conv}|^{(2)}$')
    axis padded
end

% ============================================================
% Optional: Error norm (VERY useful)
% ============================================================
figure(99); clf; set(gcf,'WindowStyle','docked');

err = x1 - x2;
err_norm = vecnorm(err,2,1);

plot(t,err_norm,'LineWidth',1.5)
grid on
title('State Error Norm')
xlabel('Time (s)')
ylabel('||x_1 - x_2||_2')


% Set all figures to light background (white) and adjust axes/text colors for visibility
figs = findall(0,'Type','figure');
for k = 1:numel(figs)
    fig = figs(k);
    % theme(fig,'light')
    theme(fig,'dark')
end

%%


% ============================================================
% Fig all states: compare each state in x1 vs x2
% Combined: three state plots per figure
% ============================================================
nStates = size(x1,1);
plotsPerFig = 3;
nFigs = ceil(nStates / plotsPerFig);

for fi = 1:nFigs
    figid = 100 + fi;
    figure(figid); clf; set(gcf,'WindowStyle','docked');
    for p = 1:plotsPerFig
        si = (fi-1)*plotsPerFig + p;
        if si > nStates, break; end
        subplot(plotsPerFig,1,p)
        plot(t, x1(si,:), '-', t, x2(si,:), '--', 'LineWidth', 1.2)
        grid on
        title(sprintf('State %d Comparison', si))
        xlabel('Time (s)')
        ylabel(sprintf('x_{%d} (units)', si))
        legend('x_1','x_2')
        axis padded
    end
    % apply theme consistently
    try
        theme(gcf,'dark')
    catch
        % ignore if theme not available
    end
end

toc