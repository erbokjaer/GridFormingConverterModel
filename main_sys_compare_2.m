% ============================================================
% Setup
% ============================================================
run setup.m
export = 0;

% run liniarize.m

run grid_and_filter_parameters.m


p = struct;

% grid/filter parameters automatically available from scripts
vars = who;
for k = 1:length(vars)
    p.(vars{k}) = eval(vars{k});
end


use_last_end_x_as_init = 1;


init_zero = 0;
init_at_x0_conditions = 0;
find_steady_state = 0;

if init_at_x0_conditions || init_zero
    save_last_end_x_as_init = 1;
else 
    save_last_end_x_as_init = 0;
end

sys_num_1 = 6;
sys_num_2 = 0;

if sys_num_1 == 0
    sys_ode_1 = @system_ode_6;
    reconstruct_fun_1 = @reconstruct_logs_6;   
elseif sys_num_1 == 1
    sys_ode_1 = @system_ode;
    reconstruct_fun_1 = @reconstruct_logs;    
elseif sys_num_1 == 2
    sys_ode_1 = @system_ode_split;
    reconstruct_fun_1 = @reconstruct_logs_split;
elseif sys_num_1 == 3
    sys_ode_1 = @system_ode_3;
    reconstruct_fun_1 = @reconstruct_logs_3;
elseif sys_num_1 == 4
    sys_ode_1 = @system_ode_4;
    reconstruct_fun_1 = @reconstruct_logs_4;
elseif sys_num_1 == 5
    sys_ode_1 = @system_ode_5;
    reconstruct_fun_1 = @reconstruct_logs_5;
    idx = setdiff(1:23, 7);
elseif sys_num_1 == 6
    sys_ode_1 = @system_ode_6;
    reconstruct_fun_1 = @reconstruct_logs_6;
elseif sys_num_1 == 7
    sys_ode_1 = @system_ode_lin;
    reconstruct_fun_1 = @reconstruct_logs_lin;
elseif sys_num_1 == 8
    sys_ode_1 = @system_ode_lin;
    reconstruct_fun_1 = @reconstruct_logs_lin_5;
elseif sys_num_1 == 9
    sys_ode_1 = @system_ode_9;
    reconstruct_fun_1 = @reconstruct_logs_9;
end

if sys_num_2 == 0
    sys_ode_2 = @system_ode_6;
    reconstruct_fun_2 = @reconstruct_logs_6;   
elseif sys_num_2 == 1
    sys_ode_2 = @system_ode;
    reconstruct_fun_2 = @reconstruct_logs;    
elseif sys_num_2 == 2
    sys_ode_2 = @system_ode_split;
    reconstruct_fun_2 = @reconstruct_logs_split;
elseif sys_num_2 == 3
    sys_ode_2 = @system_ode_3;
    reconstruct_fun_2 = @reconstruct_logs_3;
elseif sys_num_2 == 4
    sys_ode_2 = @system_ode_4;
    reconstruct_fun_2 = @reconstruct_logs_4;
elseif sys_num_2 == 5
    sys_ode_2 = @system_ode_5;
    reconstruct_fun_2 = @reconstruct_logs_5;
elseif sys_num_2 == 6
    sys_ode_2 = @system_ode_6;
    reconstruct_fun_2 = @reconstruct_logs_6;
elseif sys_num_2 == 7
    sys_ode_2 = @system_ode_lin;
    reconstruct_fun_2 = @reconstruct_logs_lin;
elseif sys_num_2 == 8
    sys_ode_2 = @system_ode_lin;
    reconstruct_fun_2 = @reconstruct_logs_lin_5;
elseif sys_num_2 == 9
    sys_ode_2 = @system_ode_9;
    reconstruct_fun_2 = @reconstruct_logs_9;
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

update_plots = 1*ones(1,100);
plotPureStates = 1;

% ============================================================
% Simulation horizon
% ============================================================
t_end = 1000;
tspan = [0 t_end];

% ============================================================
% References
% ============================================================
% Grid voltage 
Vg.y_prefault  = 1.0;    % normal voltage
Vg.y_fault     = 1.0;    % voltage during fault
Vg.y_postfault = 1.0;    % voltage after clearing

Vg.t_fault     = 400;      % fault start
Vg.t_apply     = 0.;   % ramp duration into fault

Vg.t_clear     = 401;   % clearing time
Vg.t_recover   = 0.;   % ramp duration back

% Grid phase 
Phg.t_start = 150;
Phg.t_dur   = 0;
Phg.y0      = 0;     % degrees
Phg.y1      = 0;     % degrees

% Voltage reference
V.t_start = 40;
V.t_dur   = 1;
V.y0      = 1.00;
V.y1      = 1.00;

% Reactive power reference
Q.t_start = 200;
Q.t_dur   = 0;
Q.y0      = -0.0;
Q.y1      = -0.0;


% Active power reference
P.t_start = 10;
P.t_dur   = 0;
P.y0      = 0.5;
P.y1      = 0.5;

% Mechanical power
Pm.t_start = 100;
Pm.t_dur   = 0;
Pm.y0      = 0.8;
Pm.y1      = 0.9;

vg_mag = @(t) fault_profile(t,Vg);
vg_phase_rad = @(t) ramp_signal(t, Phg.t_start, Phg.t_dur, deg2rad(Phg.y0), deg2rad(Phg.y1));
V_ref = @(t) ramp_signal(t, V.t_start,  V.t_dur,  V.y0,  V.y1);
Q_ref = @(t) ramp_signal(t, Q.t_start,  Q.t_dur,  Q.y0,  Q.y1);
P_ref     = @(t) ramp_signal(t, P.t_start,  P.t_dur,  P.y0,  P.y1);
Pm_cont   = @(t) ramp_signal(t, Pm.t_start, Pm.t_dur, Pm.y0, Pm.y1);

if init_at_x0_conditions
    vg_mag = @(t) Vg.y_prefault;
    vg_phase_rad = @(t) Phg.y0 ;
    V_ref = @(t) V.y0;
    Q_ref = @(t) Q.y0;
    P_ref = @(t) P.y0;
    Pm_cont = @(t) Pm.y0;
elseif init_zero
    vg_mag = @(t) 1;
    vg_phase_rad = @(t) 0;
    V_ref = @(t) 1;
    Q_ref = @(t) 0;
    P_ref = @(t) 0;
    Pm_cont = @(t) 0;
end

u0 = [V_ref(0); P_ref(0); Q_ref(0)];
d0 = [Pm_cont(0); vg_mag(0); vg_phase_rad(0)];
p.u0 = u0;
p.d0 = d0;

% ============================================================
% Faults and limits
% ============================================================
I_max = 1.2;

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

run temp_gains.m

% ============================================================
% Load Grid and Filter Parameters
% ============================================================

run split_system_parameters_5.m

run split_system_parameters.m




% ============================================================
% Initial state
% ============================================================

state_idx = struct( ...
        'i1d', 1, 'i1q', 2, ...         % Converter currents
        'i2d', 3, 'i2q', 4, ...         % Grid currents
        'vcd', 5, 'vcq', 6, ...         % Capacitor voltage
        'omega_g', 7, ...               % Grid speed (was 8)
        'it1d', 8, 'it1q', 9, ...       % Trap filter 1
        'vt1d', 10, 'vt1q', 11, ...
        'it2d', 12, 'it2q', 13, ...     % Trap filter 2
        'vt2d', 14, 'vt2q', 15, ...
        'delta_c', 16, 'omega_c', 17,...% Converter angle & speed
        'xi_Q', 18, ...                 % Reactive power integrator
        'xi_vd', 19, 'xi_vq', 20, ...   % Virtual impedance integrators
        'xi_id', 21, 'xi_iq', 22, ...    % Current controller integrators
        'delta_g', 23 ...
    );

if (sys_num_1 == 4) || (sys_num_1 == 5)  || (sys_num_1 == 6) || (sys_num_1 == 8) || (sys_num_1 == 9)
    x0_1 = zeros(22,1); 
else 
    x0_1 = zeros(23,1); 
end

x0_1(state_idx.vcd) = 1;
x0_1(state_idx.omega_g) = w_nom;
x0_1(state_idx.omega_c) = w_nom;
x0_1(state_idx.vt1d) = 1;
x0_1(state_idx.vt2d) = 1;

if (sys_num_2 == 4) || (sys_num_2 == 5)  || (sys_num_2 == 6) || (sys_num_2 == 8) || (sys_num_2 == 9)
    x0_2 = zeros(22,1); 
else 
    x0_2 = zeros(23,1); 
end
x0_2(state_idx.vcd) = 1;
x0_2(state_idx.omega_g) = w_nom;
x0_2(state_idx.omega_c) = w_nom;
x0_2(state_idx.vt1d) = 1;
x0_2(state_idx.vt2d) = 1;







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
checkpoint_file_1 = "init_data/system_ode_states_" + num2str(sys_num_1) + ".mat";
checkpoint_file_2 = "init_data/system_ode_states_" + num2str(sys_num_2) + ".mat";


if use_last_end_x_as_init && exist(checkpoint_file_1,'file') && not((sys_num_1 == 7) || (sys_num_1 == 8))
    load(checkpoint_file_1);
    x0_1 = last_x;
    res = norm(sys_ode_1(0, x0_1, p));
elseif (sys_num_1 == 5)
    load("init_data/system_ode_states_3.mat");
    x0_1 = last_x;
    res = norm(system_ode_3(0, x0_1, p));
elseif (sys_num_1 == 8)
    load("init_data/system_ode_states_5.mat");
    x0_1 = last_x;
    res = norm(system_ode_5(0, x0_1, p));
end
if use_last_end_x_as_init && exist(checkpoint_file_2,'file') && not((sys_num_2 == 7) || (sys_num_2 == 8))
    load(checkpoint_file_2);
    x0_2 = last_x;
    res = norm(sys_ode_2(0, x0_2, p));
elseif (sys_num_2 == 7)
    load("init_data/system_ode_states_3.mat");
    x0_2 = last_x;
    res = norm(system_ode_3(0, x0_2, p));
elseif (sys_num_2 == 8)
    load("init_data/system_ode_states_5.mat");
    x0_2 = last_x;
    res = norm(system_ode_5(0, x0_2, p));
end



% ============================================================
% Steady-state for system 1
% ============================================================

opts = optimoptions('lsqnonlin', ...
        'Display','off', ...           
        'FunctionTolerance',1e-12, ...
        'StepTolerance',1e-12, ...
        'OptimalityTolerance',1e-12, ...
        'MaxFunctionEvaluations',1e6);



if (sys_num_1 == 7)
    res1 = norm(system_ode_3(0, x0_1, p));
elseif (sys_num_1 == 8)
    res1 = norm(system_ode_5(0, x0_1, p));
else
    res1 = norm(sys_ode_1(0, x0_1, p));
end



if find_steady_state && (res1 > 1e-7)
    if (sys_num_1 == 7)
        f_ss1 = @(x) system_ode_3(0, x, p);
        x0_1 = lsqnonlin(f_ss1, x0_1, [], [], opts);
    elseif (sys_num_1 == 8)
        f_ss1 = @(x) system_ode_5(0, x, p);
        x0_1 = lsqnonlin(f_ss1, x0_1, [], [], opts);
    else
        f_ss1 = @(x) sys_ode_1(0, x, p);
        x0_1 = lsqnonlin(f_ss1, x0_1, [], [], opts);

        res1 = norm(sys_ode_1(0, x0_1, p));
        if res1 > 10000

            x0_1(18:20) = 0*x0_1(18:20);
            x0_1 = lsqnonlin(f_ss1, x0_1, [], [], opts);
        end
    end
end

% ============================================================
% Steady-state for system 2
% ============================================================
if (sys_num_2 == 7)
    res2 = norm(system_ode_3(0, x0_2, p));
elseif (sys_num_2 == 8)
    res2 = norm(system_ode_5(0, x0_2, p));
else
    res2 = norm(sys_ode_2(0, x0_2, p));
end

if find_steady_state && (res2 > 1e-7)
    if (sys_num_2 == 7)
        f_ss2 = @(x) system_ode_3(0, x, p);
        x0_2 = lsqnonlin(f_ss2, x0_2, [], [], opts);
    elseif (sys_num_2 == 8)
        f_ss2 = @(x) system_ode_5(0, x, p);
        x0_2 = lsqnonlin(f_ss2, x0_2, [], [], opts);
    else
        f_ss2 = @(x) sys_ode_2(0, x, p);
        x0_2 = lsqnonlin(f_ss2, x0_2, [], [], opts);

        res2 = norm(sys_ode_2(0, x0_2, p));
        if res2 > 10000
            x0_2(18:20) = 0*x0_2(18:20);
            x0_2 = lsqnonlin(f_ss2, 0*x0_2, [], [], opts);
        end

    end
end

% x0(state_idx.delta_c) = mod(x0(state_idx.delta_c), 2*pi) - 2*pi;
% % x0(state_idx.delta_c) = 0;
% x0(state_idx.delta_g) = mod(x0(state_idx.delta_g), 2*pi) - 2*pi;
if not((sys_num_1 == 0) || (sys_num_1 == 4) || (sys_num_1 == 5) || (sys_num_1 == 6) ) && not(sys_num_1 == 8) && not(sys_num_1 == 9)
    x0_1(state_idx.delta_c) = -mod(x0_1(state_idx.delta_g), 2*pi) + mod(x0_1(state_idx.delta_c), 2*pi);
    x0_1(state_idx.delta_g) = 0;
end

if not((sys_num_2 == 0) || (sys_num_2 == 4)  || (sys_num_2 == 5) || (sys_num_2 == 6) ) && not(sys_num_2 == 8) && not(sys_num_2 == 9)
    x0_2(state_idx.delta_c) = -mod(x0_2(state_idx.delta_g), 2*pi) + mod(x0_2(state_idx.delta_c), 2*pi);
    x0_2(state_idx.delta_g) = 0;
end




p_num = [Rlp; Kpq; Kvq; Kiq; Rv; Xv; Kpv; Kpc; L1; L2; Jconv; Jg]';





if sys_num_1 == 8
    Afull_5 = A_5 + A_fun_5(x0_1, u0, d0, p_num);
    Bfull_5 = B_5 + B_fun_5(x0_1,u0,d0,p_num);
    Efull_5 = E_5 + E_fun_5(x0_1,u0,d0,p_num);
    p.Afull = Afull_5;
    p.Bfull = Bfull_5;
    p.Efull = Efull_5;
else
    Afull = A + A_fun(x0_1, u0, d0, p_num);
    Bfull = B + B_fun(x0_1,u0,d0,p_num);
    Efull = E + E_fun(x0_1,u0,d0,p_num);
    p.Afull = Afull;
    p.Bfull = Bfull;
    p.Efull = Efull;
end




% eig(A + A_fun_5(x0, u0, d0, p_num))

p_vec = struct();

params = {'Rlp','Kpq','Kvq','Kiq','Rv','Xv','Kpv','Kpc',...
          'L1','L2','Jconv','Jg'};

for k = 1:length(params)
    p_vec.(params{k}) = sym(params{k}, 'real');
end







% ============================================================
% Solver
% ============================================================
% opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-2);
opts = odeset('RelTol',1e-4, 'AbsTol',1e-8, 'MaxStep',1e-1);

% [t,x] = ode23t(@(t,x) sys_ode(t,x,p), tspan, x0, opts);
% [t,x] = ode15s(@(t,x) sys_ode(t,x,p), tspan, x0, opts);
% [t,x] = ode45(@(t,x) sys_ode(t,x,p), tspan, x0, opts);

% t_fixed = linspace(tspan(1), tspan(end), 1000);
% [t1,x1] = ode23t(@(t,x) sys_ode_1(t,x,p), t_fixed, x0, opts);
% [t2,x2] = ode23t(@(t,x) sys_ode_2(t,x,p), t_fixed, x0, opts);


p.x0 = x0_1;
if not(sys_num_1 == 0)
    [t1,x1] = ode23t(@(t,x) sys_ode_1(t,x,p), tspan, x0_1, opts);
end

if sys_num_2 == 8
    Afull_5 = A_5 + A_fun_5(x0_2, u0, d0, p_num);
    Bfull_5 = B_5 + B_fun_5(x0_2,u0,d0,p_num);
    Efull_5 = E_5 + E_fun_5(x0_2,u0,d0,p_num);
    p.Afull = Afull_5;
    p.Bfull = Bfull_5;
    p.Efull = Efull_5;
else
    Afull = A + A_fun(x0_2, u0, d0, p_num);
    Bfull = B + B_fun(x0_2,u0,d0,p_num);
    Efull = E + E_fun(x0_2,u0,d0,p_num);
    p.Afull = Afull;
    p.Bfull = Bfull;
    p.Efull = Efull;
end


p.x0 = x0_2;
if not(sys_num_2 == 0)
    [t2,x2] = ode23t(@(t,x) sys_ode_2(t,x,p), tspan, x0_2, opts);
end

if (sys_num_1 == 0) && (sys_num_2 == 0)
    warning("No system chosen")
    return
end

if not(sys_num_1 == 0)
    t1 = t1.'; 
    x1 = x1.';   
else
    t1 = t2'*0;
    x1 = x2';
end


if not(sys_num_2 == 0)
    t2 = t2.'; 
    x2 = x2.';
else
    t2 = t1*0;
    x2 = x1;
end






if (sys_num_1 == 7) || (sys_num_1 == 8)
    p.x0 = x0_1;
    p.log0 = reconstruct_logs_split(p, 0, x0_1);
end
% p.x0 = x0;
% p.u0 = u0;
% p.d0 = d0;
log1 = reconstruct_fun_1(p,t1,x1);


if (sys_num_2 == 7) || (sys_num_2 == 8)
    p.x0 = x0_2;
    p.log0 = reconstruct_logs_split(p, 0, x0_2);
end

log2 = reconstruct_fun_2(p,t2,x2);


if save_last_end_x_as_init && not(any(isnan(x1(:, end))))

    last_x = x1(:, end);
    save(checkpoint_file_1, 'last_x');
end

if save_last_end_x_as_init && not(any(isnan(x2(:, end))))

    last_x = x2(:, end);
    save(checkpoint_file_2, 'last_x');
end


angle_idx = state_idx.delta_c;
delta0 = x0_1(angle_idx,:);

delta1 = -x1(angle_idx,:);
delta2 = -x2(angle_idx,:);
% delta1 = 0;
% delta2 = 0;

if (sys_num_1 == 7) || (sys_num_1 == 8)
    i11 = dq_rotate_lin(x1(state_idx.i1d,:), x1(state_idx.i1q,:), delta1, x0_1(state_idx.i1d,:), x0_1(state_idx.i1q,:), delta0);
    i21 = dq_rotate_lin(x1(state_idx.i2d,:), x1(state_idx.i2q,:), delta1, x0_1(state_idx.i2d,:), x0_1(state_idx.i2q,:), delta0);
    V_PCC1 = dq_rotate_lin(log1.V_PCC(1,:), log1.V_PCC(2,:), delta1, p.log0.V_PCC(1,1), p.log0.V_PCC(2,1), delta0);
    Econv1 = dq_rotate_lin(log1.Econv(1,:), log1.Econv(2,:), delta1, p.log0.Econv(1,1), p.log0.Econv(2,1), delta0);
else
    i11 = dq_rotate(x1(state_idx.i1d,:), x1(state_idx.i1q,:), delta1);
    i21 = dq_rotate(x1(state_idx.i2d,:), x1(state_idx.i2q,:), delta1);
    V_PCC1 = dq_rotate(log1.V_PCC(1,:), log1.V_PCC(2,:), delta1);
    Econv1 = dq_rotate(log1.Econv(1,:), log1.Econv(2,:), delta1);
end

delta0 = x0_2(angle_idx,:);
if (sys_num_1 == 7) || (sys_num_1 == 8)
    i12 = dq_rotate_lin(x2(state_idx.i1d,:), x2(state_idx.i1q,:), delta2, x0_2(state_idx.i1d,:),x0_2(state_idx.i1q,:),delta0);
    i22 = dq_rotate_lin(x2(state_idx.i2d,:), x2(state_idx.i2q,:), delta2, x0_2(state_idx.i2d,:), x0_2(state_idx.i2q,:), delta0);
    V_PCC2 = dq_rotate_lin(log2.V_PCC(1,:), log2.V_PCC(2,:), delta2, p.log0.V_PCC(1,1), p.log0.V_PCC(2,1), delta0);
    Econv2 = dq_rotate_lin(log2.Econv(1,:), log2.Econv(2,:), delta2, p.log0.Econv(1,1), p.log0.Econv(2,1), delta0);
else 
    i12 = dq_rotate(x2(state_idx.i1d,:), x2(state_idx.i1q,:), delta2);
    i22 = dq_rotate(x2(state_idx.i2d,:), x2(state_idx.i2q,:), delta2);
    V_PCC2 = dq_rotate(log2.V_PCC(1,:), log2.V_PCC(2,:), delta2);
    Econv2 = dq_rotate(log2.Econv(1,:), log2.Econv(2,:), delta2);
end




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
    plot(t1,i11(1,:),'-', ...
         t1,i11(2,:),'-', ...
         t2,i12(1,:),'--', ...
         t2,i12(2,:),'--','LineWidth',1)
    legend('$i_{1d}$','$i_{1q}$','$i_{1d}^{(2)}$','$i_{1q}^{(2)}$')
    grid on; title('Converter Currents'); axis padded

    subplot(3,1,2)
    plot(t1,i21(1,:),'-', ...
         t1,i21(2,:),'-', ...
         t2,i22(1,:),'--', ...
         t2,i22(2,:),'--','LineWidth',1)
    legend('$i_{2d}$','$i_{2q}$','$i_{2d}^{(2)}$','$i_{2q}^{(2)}$')
    grid on; title('Grid Side Currents'); axis padded

    subplot(3,1,3)
    plot(t1,V_PCC1(1,:),'-', ...
         t1,V_PCC1(2,:),'-', ...
         t2,V_PCC2(1,:),'--', ...
         t2,V_PCC2(2,:),'--','LineWidth',1)
    legend('$v_d$','$v_q$','$v_d^{(2)}$','$v_q^{(2)}$')
    grid on; title('PCC Voltage'); axis padded
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
    plot(t1,log1.P_conv,'-', t2,log2.P_conv,'--', ...
         t1,log1.P_grid,'-', t2,log2.P_grid,'--', ...
         t1,log1.P_PCC,'-', t2,log2.P_PCC,'--','LineWidth',1)
    legend('$P_{conv}$','$P_{conv}^{(2)}$', ...
           '$P_{grid}$','$P_{grid}^{(2)}$', ...
           '$P_{PCC}$','$P_{PCC}^{(2)}$')
    grid on; title('Active Power'); axis padded

    subplot(2,1,2)
    plot(t1,log1.Q_conv,'-', t2,log2.Q_conv,'--', ...
         t1,log1.Q_grid,'-', t2,log2.Q_grid,'--', ...
         t1,log1.Q_PCC,'-', t2,log2.Q_PCC,'--','LineWidth',1)
    legend('$Q_{conv}$','$Q_{conv}^{(2)}$', ...
           '$Q_{grid}$','$Q_{grid}^{(2)}$', ...
           '$Q_{PCC}$','$Q_{PCC}^{(2)}$')
    grid on; title('Reactive Power'); axis padded
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
    % plot(t,log1.delta_g_shifted,'-', ...
    %      t,log2.delta_g_shifted,'--', ...
    %      t,log1.delta_conv,'-', ...
    %      t,log2.delta_conv,'--','LineWidth',1.2)
    plot(t1,log1.delta_g,'-', ...
         t2,log2.delta_g,'--', ...
         t1,log1.delta_conv,'-', ...
         t2,log2.delta_conv,'--','LineWidth',1.2)
    legend('$\delta_g$','$\delta_g^{(2)}$', ...
           '$\delta_{conv}$','$\delta_{conv}^{(2)}$')
    grid on; axis padded

    subplot(2,1,2)
    plot(t1,50*log1.omega_g,'-', ...
         t2,50*log2.omega_g,'--', ...
         t1,50*log1.omega_conv,'-', ...
         t2,50*log2.omega_conv,'--','LineWidth',1.2)
    legend('$\omega_g$','$\omega_g^{(2)}$', ...
           '$\omega_{conv}$','$\omega_{conv}^{(2)}$')
    grid on; axis padded
    title('Synchronization States'); xlabel('Time (s)')
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

    plot(t1,log1.Vrefd,'-', ...
         t2,log2.Vrefd,'--', ...
         t1,log1.V_PCC_mag,'-', ...
         t2,log2.V_PCC_mag,'--','LineWidth',1.2)
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
    fig5 = findobj('Type','figure','Number',5);
    if isempty(fig5)
        fig5 = figure(5);
        set(fig5,'WindowStyle','docked');
    else
        clf(fig5);
    end
    set(0,'CurrentFigure',fig5);

    plot(t1,log1.Iref_mag_lim,'-', ...
         t2,log2.Iref_mag_lim,'--','LineWidth',1.2)
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
    fig6 = findobj('Type','figure','Number',6);
    if isempty(fig6)
        fig6 = figure(6);
        set(fig6,'WindowStyle','docked');
    else
        clf(fig6);
    end
    set(0,'CurrentFigure',fig6);

    plot(t1,log1.i1_mag,'-', t2,log2.i1_mag,'--', ...
         t1,log1.i2_mag,'-', t2,log2.i2_mag,'--','LineWidth',1.2)
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
    fig7 = findobj('Type','figure','Number',7);
    if isempty(fig7)
        fig7 = figure(7);
        set(fig7,'WindowStyle','docked');
    else
        clf(fig7);
    end
    set(0,'CurrentFigure',fig7);

    subplot(2,1,1)
    plot(t1,Econv1(1,:),'-', ...
         t1,Econv1(2,:),'-', ...
         t2,Econv2(1,:),'--', ...
         t2,Econv2(2,:),'--','LineWidth',1.2)
    legend('$E_{conv,d}$','$E_{conv,q}$', ...
           '$E_{conv,d}^{(2)}$','$E_{conv,q}^{(2)}$')
    grid on; title('Converter Internal Voltage (d/q)')
    xlabel('Time (s)'); axis padded

    subplot(2,1,2)
    plot(t1,log1.Econv_mag,'-', ...
         t2,log2.Econv_mag,'--','LineWidth',1.2)
    grid on
    title('Converter Internal Voltage Magnitude')
    xlabel('Time (s)'); ylabel('Voltage (pu)')
    legend('$|E_{conv}|$','$|E_{conv}|^{(2)}$')
    axis padded
end





fig8 = findobj('Type','figure','Number',8);
if isempty(fig8)
    fig8 = figure(8);
    set(fig8,'WindowStyle','docked');
else
    clf(fig8);
end
set(0,'CurrentFigure',fig8);

if sys_num_1 == 5 || sys_num_2 == 8
    d1 = (log1.delta_conv);
    d2 = (log2.delta_conv);
else
    d1 = (log1.delta_conv) - (log1.delta_g);
    d2 = (log2.delta_conv) - (log2.delta_g);
end

plot(t1, d1, '-', t2, d2, '--', 'LineWidth', 1.2)
grid on
title('Angle Difference: \delta_g - \delta_{conv} (x1 vs x2)')
xlabel('Time (s)')
ylabel('Angle Difference (rad)')
legend('\delta_g - \delta_{conv}','\delta_g^{(2)} - \delta_{conv}^{(2)}')
axis padded



if update_plots(9) && (sys_num_1 == 6)
    fig9 = findobj('Type','figure','Number',9);
    if isempty(fig9)
        fig9 = figure(9);
        set(fig9,'WindowStyle','docked');
    else
        clf(fig9);
    end
    set(0,'CurrentFigure',fig9);

    subplot(3,1,1)
    plot(t1,log1.V_meas,'-', t1,log1.V_meas_c,'--','LineWidth',1.2);
    grid on
    ylabel('Voltage (pu)')
    title('Measured Voltage')
    legend('V_{meas}','V_{meas}^{(c)}','Location','best')
    xlabel('Time (s)')

    subplot(3,1,2)
    plot(t1,log1.P_meas_c,'-', t1,log1.P_meas,'--','LineWidth',1.2);
    grid on
    ylabel('Active Power (pu)')
    title('Measured Active Power')
    legend('P_{meas}^{(c)}','P_{meas}','Location','best')
    xlabel('Time (s)')

    subplot(3,1,3)
    plot(t1,log1.Q_meas_c,'-', t1,log1.Q_meas,'--','LineWidth',1.2);
    grid on
    ylabel('Reactive Power (pu)')
    title('Measured Reactive Power')
    legend('Q_{meas}^{(c)}','Q_{meas}','Location','best')
    xlabel('Time (s)')


end


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

% ============================================================
% Fig all states: compare each state in x1 vs x2
% Combined: three state plots per figure
% ============================================================


if plotPureStates
    nStates1 = size(x1,1);
    nStates2 = size(x2,1);
    nStates = min(nStates1,nStates2);
    plotsPerFig = 3;
    nFigs = ceil(nStates / plotsPerFig);
    
    for fi = 1:nFigs

        figid = findobj('Type','figure','Number',100 + fi);

        if isempty(figid)
            figid = figure(100 + fi);
            set(figid,'WindowStyle','docked');
        else
            clf(figid);
        end
        set(0,'CurrentFigure',figid);


        % figure(figid); clf; set(gcf,'WindowStyle','docked');
        for p_idx = 1:plotsPerFig
            si = (fi-1)*plotsPerFig + p_idx;
            if si > nStates, break; end
            subplot(plotsPerFig,1,p_idx)
            plot(t1, x1(si,:), '-', t2, x2(si,:), '--', 'LineWidth', 1.2)
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
end


toc