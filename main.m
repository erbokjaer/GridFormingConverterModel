% ============================================================
% Setup
% ============================================================
run setup.m

use_last_end_x_as_init = 1;
save_last_end_x_as_init = 0;

init_zero = 0;
find_steady_state = 1;

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
t_end = 100;
tspan = [0 t_end];

% ============================================================
% References
% ============================================================
% Grid voltage 
Vg.y_prefault  = 1.0;    % normal voltage
Vg.y_fault     = 1;    % voltage during fault
Vg.y_postfault = 1.0;    % voltage after clearing

Vg.t_fault     = 50;      % fault start
Vg.t_apply     = 0.0;   % ramp duration into fault

Vg.t_clear     = 51;   % clearing time
Vg.t_recover   = 0.0;   % ramp duration back

% Grid phase 
Phg.t_start = 0;
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
Q.y1      = -0.;

% Active power reference
P.t_start = 10;
P.t_dur   = 0;
P.y0      = 0;
P.y1      = 1;

% Mechanical power
Pm.t_start = 10;
Pm.t_dur   = 0;
Pm.y0      = 0;
Pm.y1      = 1.1;

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
I_max = 1.2;

% ============================================================
% Swing Equation Parameters Grid
% ============================================================
J_g = 1;
D_g = 20;

% ============================================================
% Swing Equation Parameters Converter
% ============================================================
J_conv = 1;
D_conv = 70;

% ============================================================
% Reactive Power Control
% ============================================================
K_vq =  10;      
Kp_q = 10; 
tau_q_i = 0.5;
Ki_q = K_vq/tau_q_i; 

% ============================================================
% Virtual Admittance
% ============================================================
R_v = 0.15;
X_v = 0.03;
Z_v2 = sqrt(R_v^2 + X_v^2);

% Virtual impedance PI gains
Kpv = 1;
tau_v = 0.25;
Kiv = Kpv/tau_v;

% ============================================================
% Current Control
% ============================================================
Kp_c = 1; 
tau_c_i = 0.1;
Ki_c = Kp_c/tau_c_i; 

% ============================================================
% Load Grid and Filter Parameters
% ============================================================
run grid_and_filter_parameters.m


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
if find_steady_state && not(use_last_end_x_as_init)
    f_ss = @(x) system_ode(0, x, p);
    opts = optimoptions('lsqnonlin', ...
        'Display','off', ...           % can be 'none'
        'FunctionTolerance',1e-10, ...
        'StepTolerance',1e-12, ...
        'OptimalityTolerance',1e-10, ...
        'MaxFunctionEvaluations',1e5);
    try
        x_ss = lsqnonlin(@(x) f_ss(x), x0, [], [],opts);
        res = norm(f_ss(x_ss));
        if res < 1e-7
            x0 = x_ss;
        else
            warning('Steady-state solver did not converge to sufficient tolerance (residual = %g). Using original x0.', res);
        end
    catch ME
        warning(['Steady-state solver failed: ', ME.message, '. Using original x0.']);
    end
end

checkpoint_file = 'init_data/system_ode_states.mat';

if use_last_end_x_as_init && exist(checkpoint_file,'file')
    load(checkpoint_file);
    x0 = last_x;
end


% ============================================================
% Solver
% ============================================================
% opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-2);
opts = odeset('RelTol',1e-4, 'AbsTol',1e-8, 'MaxStep',1e-1);

% [t,x] = ode15s(@(t,x) system_ode(t,x,p), tspan, x0, opts);
[t,x] = ode23t(@(t,x) system_ode(t,x,p), tspan, x0, opts);
% [t,x] = ode45(@(t,x) system_ode(t,x,p), tspan, x0, opts);

x = x.';
t = t.';

% ============================================================
% Save checkpoint
% ============================================================
if save_last_end_x_as_init && not(any(isnan(x(:, end))))
    last_x = x(:, end);
    save(checkpoint_file, 'last_x');

end

% ============================================================
% Reconstruct loged signals
% ============================================================
log = reconstruct_logs(p,t,x);

% ============================================================
% Fig 1: Currents and PCC Voltage
% ============================================================
if update_plots(1)
    fig1 = findobj('Type','figure','Number',1);
    if isempty(fig1)
        fig1 = figure(1);
        set(fig1,'WindowStyle','docked');
    else
        clf(fig1);
    end
    set(0,'CurrentFigure',fig1);
    
    subplot(3,1,1)
    plot(t,x(state_idx.i1d,:),t,x(state_idx.i1q,:),'LineWidth',1)
    legend('$i_{1d}$','$i_{1q}$','Location','best')
    grid on
    title('Converter Currents')
    axis padded
    
    subplot(3,1,2)
    plot(t,x(state_idx.i2d,:),t,x(state_idx.i2q,:),'LineWidth',1)
    legend('$i_{2d}$','$i_{2q}$','Location','best')
    grid on
    title('Grid Side Currents')
    axis padded
    
    subplot(3,1,3)
    plot(t,log.V_PCC(1,:),t,log.V_PCC(2,:),'LineWidth',1)
    legend('$v_{d}$','$v_{q}$','Location','best')
    grid on
    title('PCC Voltage')
    axis padded
end

% ============================================================
% Fig 2: Active and Reactive Power
% ============================================================
if update_plots(2)
    fig2 = findobj('Type','figure','Number',2);
    if isempty(fig2)
        fig2 = figure(2);
        set(fig2,'WindowStyle','docked');
    else
        clf(fig2);
    end
    set(0,'CurrentFigure',fig2);
    
    subplot(2,1,1)
    plot(t,log.P_conv,t,log.P_grid,t,log.P_PCC,'LineWidth',1)
    legend('$P_{conv}$','$P_{grid}$','$P_{PCC}$','Location','best')
    grid on
    title('Active Power')
    axis padded
    
    subplot(2,1,2)
    plot(t,log.Q_conv,t,log.Q_grid,t,log.Q_PCC,'LineWidth',1)
    legend('$Q_{conv}$','$Q_{grid}$','$Q_{PCC}$','Location','best')
    grid on
    title('Reactive Power')
    axis padded
end

% ============================================================
% Fig 3: Synchronization States
% ============================================================
if update_plots(3)
    fig3 = findobj('Type','figure','Number',3);
    if isempty(fig3)
        fig3 = figure(3);
        set(fig3,'WindowStyle','docked');
    else
        clf(fig3);
    end
    set(0,'CurrentFigure',fig3);
    
    subplot(2,1,1)
    plot(t,log.delta_g_shifted,'LineWidth',1.2)
    hold on
    plot(t,log.delta_conv,'LineWidth',1.2)
    hold off
    legend('$\delta_g$','$\delta_{conv}$','Location','best')
    grid on
    axis padded
    
    subplot(2,1,2)
    plot(t,log.omega_g,'LineWidth',1.2)
    hold on
    plot(t,log.omega_conv,'LineWidth',1.2)
    hold off
    legend('$\omega_g$','$\omega_{conv}$','Location','best')
    grid on
    axis padded
    title('Synchronization States')
    xlabel('Time (s)')
end

% ============================================================
% Fig 4: Vref and PCC Voltage
% ============================================================
if update_plots(4)
    fig4 = findobj('Type','figure','Number',4);
    if isempty(fig4)
        fig4 = figure(4);
        set(fig4,'WindowStyle','docked');
    else
        clf(fig4);
    end
    set(0,'CurrentFigure',fig4);
    
    plot(t,log.Vrefd,'LineWidth',1.2)
    hold on
    plot(t,log.V_PCC_mag,'--','LineWidth',1)
    hold off
    grid on
    title('Voltage Reference and PCC Voltage')
    xlabel('Time (s)')
    ylabel('Voltage (pu)')
    legend('$V_{ref,d}$','$V_{PCC}$','Location','best')
    axis padded
end

% ============================================================
% Fig 5: Reference Current Magnitude
% ============================================================
if update_plots(5)
    fig5 = findobj('Type','figure','Number',5);
    if isempty(fig5)
        fig5 = figure(5);
        set(fig5,'WindowStyle','docked');
    else
        clf(fig5);
    end
    set(0,'CurrentFigure',fig5);
    
    plot(t,log.Iref_mag_lim,'LineWidth',1.2)
    grid on
    title('Reference Current Magnitude')
    xlabel('Time (s)')
    ylabel('$I_{ref}$ (pu)')
    axis padded
end

% ============================================================
% Fig 6: Current Magnitudes
% ============================================================
if update_plots(6)
    fig6 = findobj('Type','figure','Number',6);
    if isempty(fig6)
        fig6 = figure(6);
        set(fig6,'WindowStyle','docked');
    else
        clf(fig6);
    end
    set(0,'CurrentFigure',fig6);
    
    plot(t,log.i1_mag,'LineWidth',1.2)
    hold on
    plot(t,log.i2_mag,'LineWidth',1.2)
    hold off
    grid on
    title('Current Magnitudes')
    xlabel('Time (s)')
    ylabel('Current (pu)')
    legend('$|i_1|$','$|i_2|$','Location','best')
    axis padded
end

% ============================================================
% Fig 7: Converter Internal Voltage (dq and magnitude)
% ============================================================
if update_plots(7)
    fig7 = findobj('Type','figure','Number',7);
    if isempty(fig7)
        fig7 = figure(7);
        set(fig7,'WindowStyle','docked');
    else
        clf(fig7);
    end
    set(0,'CurrentFigure',fig7);
    
    subplot(2,1,1)
    plot(t, log.Econv(1,:), t, log.Econv(2,:), 'LineWidth', 1.2)
    legend('$E_{conv,d}$','$E_{conv,q}$','Location','best')
    grid on
    title('Converter Internal Voltage (d/q)')
    xlabel('Time (s)')
    axis padded
    
    subplot(2,1,2)
    plot(t, log.Econv_mag, 'LineWidth', 1.2)
    grid on
    title('Converter Internal Voltage Magnitude')
    xlabel('Time (s)')
    ylabel('Voltage (pu)')
    legend('$|E_{conv}|$','Location','best')
    axis padded
end
% ============================================================
% Export Figures
% ============================================================
fig1Name = outputFolder + "\Currens_and_PCC_Voltage_dq" + fig_type;
fig2Name = outputFolder + "\Active_and_Reactive_Power" + fig_type;
fig3Name = outputFolder + "\Synchronization_States" + fig_type;
fig4Name = outputFolder + "\Vref_and_Capacitor_Voltage" + fig_type;
fig5Name = outputFolder + "\Iref_Magnitude" + fig_type;
fig6Name = outputFolder + "\I1_and_I2_Magnitude" + fig_type;
fig7Name = outputFolder + "\Econv_and_EconvMag" + fig_type;

if export
    exportLight(fig1,fig1Name)
    exportLight(fig2,fig2Name)
    exportLight(fig3,fig3Name)
    exportLight(fig4,fig4Name)
    exportLight(fig5,fig5Name)
    exportLight(fig6,fig6Name)
    exportLight(fig7, fig7Name);
end

if export
    % update all legend locations to specific locations (individual per figure) and export
    locs = {'west','east','east','southeast','southeast','southeast'}; % set desired locations per figure
    figs = {fig1,fig2,fig3,fig4,fig5,fig6};
    names = {fig1Name,fig2Name,fig3Name,fig4Name,fig5Name,fig6Name};

    for k = 1:numel(figs)
        fig = figs{k};
        if isvalid(fig)
            axs = findall(fig,'Type','axes');
            for a = 1:numel(axs)
                lg = legend(axs(a));
                if ~isempty(lg) && isvalid(lg)
                    lg.Location = locs{k};
                end
            end
        end
    end

    % export all
    for k = 1:numel(figs)
        exportLight(figs{k}, names{k});
    end
end

toc