function res = example_nonlinearDA_reach_09_SMIB_adaptive()
% example_nonlinearDA_reach_09_SMIB_adaptive - example of nonlinear
%    differential-algebraic reachability analysis
%
% Syntax:
%    example_nonlinearDA_reach_09_SMIB_adaptive()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       30-August-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% SMIB system with 3 diff. variables and 8 alg. variables

dim_x = 3;
dim_y = 8;

[P,x0,y0] = aux_aux3_modelParameters();
params.P = P;

% Parameters --------------------------------------------------------------

params.x0 = x0;
params.y0guess = y0.';
Bound_x = [0.005 0 0; 0 0.5e-4 0; 0 0 0.005];

% initial set
params.R0 = zonotope([params.x0.',Bound_x]); % x0

% uncertain inputs
u0 = [1.2242; 0; 1.0241];
Bound_u = zeros(dim_x,1);        % Uncertainty
params.U=zonotope([u0,Bound_u]); % Initial input state


% Reachability Settings ---------------------------------------------------

options.alg = 'lin-adaptive';
options.verbose = true;
options.tensorOrder = 2;


% Reachability Analysis ---------------------------------------------------

% Normal mode
P.mode = 'normal';
dynHandle = @(x,y,u) SMIBdyn(x,y,u,P);
conHandle = @(x,y,u) SMIBcon(x,y,u,P);

sys{1} = nonlinDASys('normal',dynHandle,conHandle);

% Fault mode
P.mode = 'fault';
dynHandle = @(x,y,u) SMIBdyn(x,y,u,P);
conHandle = @(x,y,u) SMIBcon(x,y,u,P);

sys{2} = nonlinDASys('fault',dynHandle,conHandle);
    

tswitch = [0 0.1 0.3 1];
modes = {'normal','fault','normal'};
syshandles = {sys{1},sys{2},sys{1}};

tic;
% loop over all mode switches
for i = 1:length(modes)

    % update time
    params.tStart = tswitch(i);
    params.tFinal = tswitch(i+1);
    
    % compute reachable set
    Rtemp = reach(syshandles{i}, params, options);
   
    if i == 1
    	R = Rtemp;  
    else
    	R = add(R,Rtemp); 
    end
    
    % Initial state of last state of previous mode
    params.R0 = Rtemp.timePoint.set{end};
end
compTime = toc;


% Simulation --------------------------------------------------------------

tswitch = [0 0.1 0.3 1];
modes = {'normal','fault','normal'};
syshandles = {sys{1},sys{2},sys{1}};

runs = 30;
Bound_x = [0.005 0 0; 0 0.5e-4 0; 0 0 0.005];
params.R0 = zonotope([params.x0.',Bound_x]);
dim_x = dim(params.R0);

t = cell(runs,1); x = cell(runs,1);
for r=1:runs
    % set init conditions for run
    if r<=20
        params.u = randPoint(params.U,1,'extreme');
        params.x0 = randPoint(params.R0,1,'extreme');
    else
        params.u = randPoint(params.U);
        params.x0 = randPoint(params.R0);
    end
    for i=1:length(modes)
%         P.mode = modes{i};
        params.tStart = tswitch(i);
        params.tFinal = tswitch(i+1);
        [t_new,x_new] = simulate(syshandles{i},params);
        t{r} = [t{r}; t_new]; x{r} = [x{r}; x_new];
        % update for new mode
        xFinal = x{r}(end,:)';
        params.x0 = xFinal(1:dim_x);
    end
end

% write to simResult object
simRes = simResult(x,t);


% Visualization -----------------------------------------------------------

% post-processing of R for plots
theta = {};
P = params.P;

tVec = 0;

for j = 1:size(R,1)
    for k = 1:length(R(j,1).timeInterval.set)
        
        Rcont = R(j,1).timeInterval.set{k};
        Ry = R(j,1).timeInterval.algebraic{k};
        
        Xhull      = interval(Rcont);
        delta_para = interval(infimum(Xhull(1)),supremum(Xhull(1)));
        Yhull      = interval(Ry);
        id         = interval(infimum(Yhull(3)),supremum(Yhull(3)));
        iq         = interval(infimum(Yhull(4)),supremum(Yhull(4)));

        theta{end+1,1} = (P.Tm-(P.xq-P.x_d)*id*iq) / (2*P.H*delta_para);
        theta{end,2} = iq;
        theta{end,3} = id/delta_para;
        
        tVec(end+1,1) = supremum(R(j,1).timeInterval.time{k});
    end
end


% first plot
figure;
for p=1:size(theta,2)
    subplot(3,1,p); hold on; box on;
    useCORAcolors("CORA:contDynamics");
    color = CORAcolor('CORA:next');
    curr_theta = theta(:,p);
    for i=1:length(curr_theta)
        IH = interval([infimum(curr_theta{i}); infimum(curr_theta{i})],...
            [supremum(curr_theta{i});supremum(curr_theta{i})]);
        IH = interval([tVec(i);infimum(IH(1))], [tVec(i+1);supremum(IH(1))]);
        plot(IH,[1 2],'FaceColor',color);
    end 
    % label
    xlabel('time')
    ylabel(p)
end

% second plot
figure; hold on; box on;

projDim = [1, 2];

% plot reachable sets
useCORAcolors("CORA:contDynamics")
plot(R,projDim,'DisplayName','Reachable set');

% plot initial set
plot(R(1).R0,projDim, 'DisplayName','Initial set');

% plot simulation results      
plot(simRes,projDim,'DisplayName','Simulations');

legend()


% examples completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function [P,x0,y0] = aux_aux3_modelParameters()
% define parameter of the model

    % Assign fixed parameters ---------------------------------------------
    P.H         = 3.5;   	% Intertia coeff.
    P.D         = 0;        % Damping constant
    P.T_d0      = 8;        % Time constant d_axis
    P.T_qO      = 1;        % Time constant q_axis
    P.xd        = 1.81;  	% d-Achsen synchrone Reaktanz
    P.x_d       = 0.3;  	% d-Achsen transiente Reaktanz  
    P.xq        = 1.76;
    P.x_q       = 0.65;
    P.ra        = 0.003;
    P.omegaS    = 2*pi*50;

    % Network parameters --------------------------------------------------
    P.mode      = 'normal';    % execution mode
    P.xT        = 0.15;        % Reactance of the transformers
    P.xl1       = 0.5;         % Reactance of the 1. line
    P.xl2       = 0.93;        % Reactance of the 2. line

    % Resulting overall network reactance for normal and fault case -------
    P.xs        = P.xT + 1/(1/P.xl1 + 1/P.xl2); % Pre-fault
    P.xs2       = P.xT; % Fault
    P.xs3       = P.xT + P.xl1; % Post-fault

    % Power Flow: Initial values ------------------------------------------
    % Bus-1 is PV-bus
    P.v1        = 1;           % voltage at bus-1
    P.p1        = 0.9;         % active power at bus-1

    % Bus-2 is slack bus
    P.v2        = 0.90081;     % Voltage at Infinite Bus := const
    P.theta2    = 0;           % Phase at Infinite Bus := const

    % Set options for nonlinear solver to obtain initial variables of
    % remaining variables of the system
    options_fsolve = optimoptions('fsolve','Display','off');
    X0netz = fsolve(@(X)aux_aux3_init_sim_netzwerk(X,P),ones(4,1),options_fsolve);

    P.theta1  = X0netz(1);     % phase at bus-1
    P.q1      = X0netz(2);     % reactive power at bus-1

    % Generator Initinal values -------------------------------------------
    X0gen=fsolve(@(X)aux_aux3_init_generator(X,P),ones(10,1),options_fsolve);

    % diff. variables
    delta0  = X0gen(1);        % state variable
    omega0  = X0gen(2);        % state variable
    Eq_0    = X0gen(3);        % state variable

    x0 = [delta0 omega0 Eq_0];

    % alg. variables
    ed0     = X0gen(4);        % alg. Variable
    eq0     = X0gen(5);        % alg. Variable
    Ef      = X0gen(6);        % Constant field voltage
    id0     = X0gen(7);        % alg. Variable
    iq0     = X0gen(8);        % alg. Variable
    P.Tm    = X0gen(9);        % Constant mechanical torque
    P.Te    = X0gen(10);       % alg. Variable

    y0 = [ed0 eq0 id0 iq0 P.p1 P.q1 P.theta1 P.v1];
end

function X0 = aux_aux3_init_sim_netzwerk(X,P)
% Initialisation of the grid, POWER FLOW
    theta1 = X(1);
    q1     = X(2);
    p2     = X(3);
    q2     = X(4);

    X0=[P.v1*P.v2/P.xs*sin(theta1-P.theta2)-P.p1;...       %Netzwerk
        P.v1*P.v2/P.xs*sin(P.theta2-theta1)-p2;...
        P.v1*P.v1/P.xs-P.v1*P.v2/P.xs*cos(theta1-P.theta2)-q1;...
        P.v2*P.v2/P.xs-P.v1*P.v2/P.xs*cos(P.theta2-theta1)-q2;...
    ];
end

function X0 = aux_aux3_init_generator(X,P)
% Assignment of the variables

    delta=X(1);
    omega=X(2);
    Eq_=X(3);
    ed=X(4);
    eq=X(5);
    Ef=X(6);
    id=X(7);
    iq=X(8);
    Tm=X(9);
    Te=X(10);

    X0=[P.omegaS*omega;...            
    (Tm-Te)/(2*P.H);...
    (Ef-Eq_-(P.xd-P.x_d)*id)/P.T_d0;...
    -ed+P.xq*iq-P.ra*id;...          
    -eq+Eq_-P.ra*iq-P.x_d*id;...
    -P.p1+ed*id+eq*iq;...       
    -P.q1+eq*id-ed*iq;...
    -Te+(eq+P.ra*iq)*iq+(ed+P.ra*id)*id;...     
    -ed+P.v1*sin(delta-P.theta1);...  
    -eq+P.v1*cos(delta-P.theta1);...
    ];
end

% ------------------------------ END OF CODE ------------------------------
