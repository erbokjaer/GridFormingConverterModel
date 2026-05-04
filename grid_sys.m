% ============================================================
% Setup
% ============================================================
run setup.m
export = 0;

% run liniarize.m

run grid_and_filter_parameters.m


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
P.y0      = 0.0;
P.y1      = 0.5;

% Mechanical power
Pm.t_start = 100;
Pm.t_dur   = 0;
Pm.y0      = 0.0;
Pm.y1      = 0.5;

vg_mag = @(t) fault_profile(t,Vg);
vg_phase_rad = @(t) ramp_signal(t, Phg.t_start, Phg.t_dur, deg2rad(Phg.y0), deg2rad(Phg.y1));
V_ref = @(t) ramp_signal(t, V.t_start,  V.t_dur,  V.y0,  V.y1);
Q_ref = @(t) ramp_signal(t, Q.t_start,  Q.t_dur,  Q.y0,  Q.y1);
P_ref     = @(t) ramp_signal(t, P.t_start,  P.t_dur,  P.y0,  P.y1);
Pm_cont   = @(t) ramp_signal(t, Pm.t_start, Pm.t_dur, Pm.y0, Pm.y1);


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


A = [...
-(R1+Rlp)/L1, 0,              Rlp/L1,        0,             -1/L1, 0,              Rlp/L1,        0, 0, 0,  Rlp/L1,        0, 0, 0, 0;
0,        -(R1+Rlp)/L1, 0,             Rlp/L1,        0,     -1/L1,          0,             Rlp/L1, 0, 0,  0,             Rlp/L1, 0, 0, 0;

Rlp/L2,        0,              -(Rlp+R2)/L2, 0,              1/L2,  0,             -Rlp/L2,       0, 0, 0, -Rlp/L2,       0, 0, 0, 0;
0,              Rlp/L2,        0,             -(Rlp+R2)/L2,   0,     1/L2,          0,             -Rlp/L2, 0, 0,  0,             -Rlp/L2, 0, 0, 0;

1/Clp,         0,             -1/Clp,        0,              0,     0,             -1/Clp,        0, 0, 0, -1/Clp,        0, 0, 0, 0;
0,              1/Clp,        0,             -1/Clp,         0,     0,              0,            -1/Clp, 0, 0,  0,            -1/Clp, 0, 0, 0;

Rlp/Lt1,       0,             -Rlp/Lt1,      0,              1/Lt1, 0,             -(Rlp+Rt1)/Lt1,0, -1/Lt1, 0, -Rlp/Lt1, 0, 0, 0, 0;
0,              Rlp/Lt1,       0,            -Rlp/Lt1,       0,     1/Lt1,          0,            -(Rlp+Rt1)/Lt1, 0, -1/Lt1, 0, -Rlp/Lt1, 0, 0, 0;

0,              0,             0,             0,              0,     0,              1/Ct1,        0, 0, 0, 0, 0, 0, 0, 0;
0,              0,             0,             0,              0,     0,              0,            1/Ct1, 0, 0, 0, 0, 0, 0, 0;

Rlp/Lt2,       0,             -Rlp/Lt2,      0,              1/Lt2, 0,             -Rlp/Lt2,      0, 0, 0, -(Rlp+Rt2)/Lt2, 0, -1/Lt2, 0, 0;
0,              Rlp/Lt2,       0,            -Rlp/Lt2,       0,     1/Lt2,          0,            -Rlp/Lt2, 0, 0, 0, -(Rlp+Rt2)/Lt2, 0, -1/Lt2, 0;

0,              0,             0,             0,              0,     0,              0,            0, 0, 0, 1/Ct2, 0, 0, 0, 0;
0,              0,             0,             0,              0,     0,              0,            0, 0, 0, 0, 1/Ct2, 0, 0, 0;
0,              0,             0,             0,              0,     0,              0,            0, 0, 0, 0, 0, 0, 0, -Dg/Jg];
% e1 = (eig(A))
% 
% A(15,3) = 1/Jg;
% 
% e2 = (eig(A))

T = [0, w_nom;-w_nom,0];
J = blkdiag(T,T,T,T,T,T,T,0);
Agem = A + J;
% Agem = A;
real(eig(Agem));



Bgem = [...
1/L1, 0;
0,    1/L1;
0,      0;
0,      0;
0,     0;
0,     0;
0,     0;
0,     0;
0,     0;
0,     0;
0,     0;
0,     0;
0,     0;
0,     0;
0,     0];

Egem = [...
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
0,     0, 0;
-1/Jg,     0, 0];

C_gem = [...
Rlp, 0,  -Rlp, 0, 1, 0,  -Rlp, 0, 0, 0,  -Rlp, 0, 0, 0;
0,  Rlp, 0, -Rlp, 0, 1,   0, -Rlp, 0, 0,   0, -Rlp, 0, 0;
0,  0,   1, 0,    0, 0,   0, 0,    0, 0,   0, 0,    0, 0;
0,  0,   0, 1,    0, 0,   0, 0,    0, 0,   0, 0,    0, 0];


% Add one state to all system matrices/vectors (to make state dimension 16)
% Existing matrices: A, Agem, Bgem, Egem, C_gem, J (blkdiag), and any vectors sized to 15
nOld = size(A,1);         % expected 15
nNew = nOld + 1;          % 16

if nOld ~= 15
    warning('Expected original state dimension 15, got %d. Proceeding to add one state.', nOld);
end

% Expand A: add zero row and column (new state has zero linear coupling by default)
A(nNew,nNew) = 0;
A(nNew,1:nNew-1) = 0;
A(1:nNew-1,nNew) = 0;

% Expand Agem similarly (Agem was A + J)
Agem(nNew,nNew) = 0;
Agem(nNew,1:nNew-1) = 0;
Agem(1:nNew-1,nNew) = 0;

% Expand Bgem: add zero row (no direct input to new state)
if size(Bgem,1) < nNew
    Bgem(nNew,size(Bgem,2)) = 0;
end

% Expand Egem: add zero row (no direct disturbance to new state)
if size(Egem,1) < nNew
    Egem(nNew,size(Egem,2)) = 0;
end

% Expand C_gem: add zero column for the new state if needed (C maps states to outputs)
if size(C_gem,2) < nNew
    C_gem(:,nNew) = 0;
end

% Expand J (used in fnl): J should be nNew-by-nNew skew block; pad with zeros
if size(J,1) < nNew
    J(nNew,nNew) = 0;
    J(nNew,1:nNew-1) = 0;
    J(1:nNew-1,nNew) = 0;
end

% Ensure x0, x_star size compatibility where they exist later
% If x0 exists already, pad to nNew
if exist('x0','var') && numel(x0) < nNew
    x0(nNew,1) = 0;
end

% If x_star exists, pad as well
if exist('x_star','var') && numel(x_star) < nNew
    x_star(nNew,1) = 0;
end

% Update any hard-coded sizes used elsewhere (tspan, plotting) by setting a variable
nStates = nNew;





p = struct;

% grid/filter parameters automatically available from scripts
vars = who;
for k = 1:length(vars)
    p.(vars{k}) = eval(vars{k});
end



x0 = zeros(16,1);

x0(5) = 1;
x0(9) = 1;
x0(13) = 1;
% x0(15) = 1;

% x0 = [
%     0.4541;
%     0.2085;
%     0.4578;
%     0.1896;
%     0.9874;
%     0.1902;
%    -0.0008;
%     0.0042;
%     0.9864;
%     0.1975;
%    -0.0004;
%     0.0021;
%     0.9861;
%     0.1975;
%    0;
% ];



delta = 0.;
u_star = [cos(delta); sin(delta)];


% sqrt(4*Dg*R2)

d_star = [0; 1; 0.];
x0(1) = d_star(1);
x0(3) = d_star(1);




% options = optimoptions('fsolve','Display','off','TolFun',1e-12);
% x_star = fsolve(@(x) gem_system(0, x, u_star, d_star, p), x0, options);
opts = optimoptions('lsqnonlin', ...
        'Display','off', ...           
        'FunctionTolerance',1e-12, ...
        'StepTolerance',1e-12, ...
        'OptimalityTolerance',1e-12, ...
        'MaxFunctionEvaluations',1e6);
x_star = lsqnonlin(@(x) gem_system(0, x, u_star, d_star, p), x0, [], [], opts);

x0 = x_star;

delta = 0.;
% u = 1*[cos(delta) * 1 ; sin(delta)*1];
% d = [1; 1000000/2; 0];
u = u_star;
d = d_star;

% x0 = zeros(15,1);
function dx = gnl(x,d,p)
dx = zeros(16,1);
dx(3) = -(1/p.L2) * d(2) * cos(d(3));
dx(4) = -(1/p.L2) * d(2) * sin(d(3));
dx(15) = (1/p.Jg) * (d(2)*cos(d(3))*x(3)+ d(2)*sin(d(3))*x(4) );    
end

function dx = fnl(x,p)
dx = x(16).*p.J*x;    
end
% Agem*x_star + fnl(x_star,p)
% fnl(x_star,p)
% gnl(x_star,d_star,p)
% Bgem*u_star 
% Egem*d_star
res = Agem*x_star + fnl(x_star,p) + gnl(x_star,d_star,p) + Bgem*u_star + Egem*d_star

% dx_tilde = Agem*x_tilde + fnl(x_star + x_tilde) - fnl(x_star) + gnl(x_star + x_tilde, d_star + d_tilde) - gnl(x_star, d_star) + Bgem*(u - u_star) + Egem*(d - d_star);

% dx_tilde = Agem*x_tilde + x(15).*p.J*x_tilde + gnl(x_star + x_tilde, d_star + d_tilde) - gnl(x_star, d_star) + Bgem*(u - u_star) + Egem*(d - d_star);

% x0 = zeros(16,1);
% x0(5) = 1;
% x0(9) = 1;
% x0(13) = 1;
% x0(15) = 0;


% delta = 0.25;
% u = [cos(delta) ; sin(delta)];
% d = [0; 0; 0.];



[t,x] = ode23t(@(t,x) gem_system(t,x,u,d,p), tspan, x0);



t = t';
x = x';
sqrt(x(1,end)^2 + x(1,end));
x(15,:) =x(15,:)  * 50;



% Plot states: each plot (subplot) should show two states, and each figure max three subplots
nStates = size(x,1);
statesPerPlot = 2;
maxSubplotsPerFigure = 3;
plotsPerFigure = maxSubplotsPerFigure; % each subplot contains two states

% Define state names (modify or extend if state count changes)
defaultStateNames = arrayfun(@(i) sprintf('x_{%d}', i), 1:nStates, 'UniformOutput', false);
% Provide more meaningful names when known (example mapping for 15 states)
if nStates >= 15
    stateNames = {...
        'i1_d','i1_q',...       % 1-2
        'i2_d','i2_q',...       % 3-4
        'v_lp_d','v_lp_q',...   % 5-6
        'i_t1_d','i_t1_q',...   % 7-8
        'v_t1_d','v_t1_q',...   % 9-10
        'i_t2_d','i_t2_q',...   % 11-12
        'v_t2_d','v_t2_q',...   % 13-14
        'omega_g'...            % 15
        };
    % If more states than 15, append defaults
    if numel(stateNames) < nStates
        stateNames = [stateNames, defaultStateNames(numel(stateNames)+1:nStates)];
    end
else
    stateNames = defaultStateNames;
end

% Determine total number of subplots needed
nSubplots = ceil(nStates / statesPerPlot);
nFigures = ceil(nSubplots / plotsPerFigure);

stateIdx = 1;
baseFig = 1; % starting figure number; keep figure numbers consistent across runs
for f = 1:nFigures
    figNum = baseFig + (f-1); % preserve and increment figure numbers explicitly
    figure(figNum);
    clf(figNum); % clear figure but keep number
    for sp = 1:plotsPerFigure
        subplotIdx = (f-1)*plotsPerFigure + sp;
        if subplotIdx > nSubplots
            break;
        end
        subplot(plotsPerFigure,1,sp);
        hold on;
        legendEntries = {};
        % plot up to two states in this subplot
        firstIdx = stateIdx;
        for k = 1:statesPerPlot
            if stateIdx > nStates
                break;
            end
            plot(t, x(stateIdx,:), 'LineWidth', 1.2);
            legendEntries{end+1} = stateNames{stateIdx}; %#ok<SAGROW>
            stateIdx = stateIdx + 1;
        end
        grid on;
        ylabel('Amplitude');
        title(sprintf('States %d to %d', firstIdx, stateIdx - 1));
        if sp == plotsPerFigure || subplotIdx == nSubplots
            xlabel('Time (s)');
        end
        legend(legendEntries, 'Location', 'best');
        hold off;
    end
end