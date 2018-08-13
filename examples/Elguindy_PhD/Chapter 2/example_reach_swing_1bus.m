function completed = example_reach_swing_1bus()

% Syntax:  
%    Reach_Swing_1bus
%
% Inputs:
%    no
%
% Outputs:
%    no 
%
% Example: 
% 
% 
% Author:        Ahmed El-Guindy
% Written:       11-November-2015
% Last update:   ---
% Last revision: 28-August-2017

%------------- BEGIN CODE --------------


% Define a global variable for SMIB_Swing_Const.m 
global P I

% SMIB system with 3 diff. variables and 8 alg. variables
dim_x = 2;
dim_y = 4;

options.tensorOrder = 2;
reachflag = 1;

% Set options -------------------------------------------------------------
options.tStart   = 0;      % Start time
options.timeStep = 0.0001;   % Time step size for reachable set computation

%% Assign fixed parameters ------------------------------------------------
%  Example 13.2 p. 864
P.M1         = 1/(15*pi);   	% Tr???gheitskoeffizient
P.D1         = 0.5;           % D???mpfungskoeffizient
P.xd1       = 1/(1i*0.2);  	% d-Achsen transiente Reaktanz     
P.omegaS     = 2*pi*50;

%% Network parameters
P.xT        = 1*0.15;     % Reactance of the transformer at Bus 1
P.xl1       = 1*0.5;      % Reactance of trasmission line 1
P.xl2       = 1*0.93;     % Reactance of trasmission line 2

%% Resulting overall network reactance for normal and fault case
%  Series and parallel computation  
P.xs        = P.xT + 1/(1/P.xl1 + 1/P.xl2); % No fault
P.xs3       = P.xT + P.xl1;                 % Fault in line 2
% Fault sequence (This can be changed according to the needs)

%% Initial values for the Power Flow
% Bus 1 is a PV-bus (connected to a generator)
% Bus 2 is an infinite bus
P.v1        = 1;           % Voltage at bus 1
P.p1        = 0.9;         % Active power at bus 1
P.v2        = 0.90081;     % Voltage at Infinite Bus := const
P.theta2    = 0;           % Phase at Infinite Bus   := const

%% Compute initial values for the SMIB system
options_fsolve = optimoptions('fsolve','Display','off'); % Turn off display

% Use fsolve to get theta1 and q1
% Phase angle and reactive power at bus 1 (Unknown for PV-bus)
X0net       = fsolve(@(X)init_network(X, P), ones(2,1),options_fsolve);
I.theta1    = X0net(1);
I.Q1        = X0net(2);      

I.P1        = 0.9;
I.V1        = 1;

% Use fsolve to get synchrous generator paramters
X0gen       = fsolve(@(X)init_generator(X, P,I), ones(4,1),options_fsolve);
%Zuweisung der Variablen
delta0      = X0gen(1);        % rotor angle (diff. variable)
omega0      = X0gen(2);        % rotor angular speed (diff. variable)
I.E1        = X0gen(3);        % q-axis transient voltage (diff. variable)
I.Pm1       = X0gen(4);        % d-axis machine voltage (alg. varaible)

x0 = [X0gen(1);X0gen(2)];
y0 = [I.P1 I.Q1 I.theta1 I.V1];

%% Assign parameters for reachability analysis ----------------------------

options.x0 = x0;
options.y0 = y0; 
options.y0guess = y0.';

Bound_x = 0.8e-3*eye(2);

options.R0  = zonotope([options.x0,Bound_x]); % x0
options.R_y = zonotope([options.y0.',zeros(dim_y,dim_y)]); % y0

%% Initial set for the input variables ------------------------------------
options.uTrans = 0;
Bound_u(1,1)   = 0;                       % Uncertainty
options.U=zonotope([zeros(1,1),Bound_u]); % Initial input state
options.originContained = 1;              % No origin

options.zonotopeOrder = 300; % Zonotope order
options.polytopeOrder = 10;  % Polytope order

options.taylorTerms = 10;

options.tensorOrder          = 2;
options.advancedLinErrorComp = 1;
options.intermediateOrder    = 3;

options.reductionTechnique = 'girard';

options.errorOrder = 3;
options.maxError   = 0;
options.oldError_x = 0.001*ones(dim_x,1);
options.maxError_x = options.oldError_x;
options.oldError_y = 0.005*ones(dim_y,1);
options.maxError_y = options.oldError_y;

options.reductionInterval = 1e5;
%--------------------------------------------------------------------------

% Specify continuous dynamics----------------------------------------------
%save powerDyn_SwingFault powerDyn_SwingFault

options.mode = 'normal'; % No fault
P.mode       = 'normal';
powerDyn_SwingNormal = nonlinDASys(2,4,1,@SMIB_Swing,@SMIB_Swing_Const,options);
save powerDyn_SwingNormal powerDyn_SwingNormal

%--------------------------------------------------------------------------
% Compute reachable set
if reachflag
    %% normal operation
    options.tStart = 0;   % Start time
    options.tFinal = 0.01; % Final time
    tic
    [Rcont1,~,Rcont_y1] = reach(powerDyn_SwingNormal, options);
    tComp = toc;
    disp(['computation time of reachable set: ',num2str(tComp)]);

    
    %% Fault
    options.mode = 'fault'; % No fault
    P.mode       = 'fault';
    powerDyn_SwingFault = nonlinDASys(2,4,1,@SMIB_Swing,@SMIB_Swing_Const,options);
    options.tStart = 0.01; % Start time (tFinal of normal operation)
    options.tFinal = 0.02; % Final time
    options.R0   = Rcont1{length(Rcont1)}{1}; % Initial state of last state from normal operation
    options.R0_y = Rcont_y1{length(Rcont_y1)}{1};
    [Rcont2,~,Rcont_y2] = reach(powerDyn_SwingFault, options);

    %% normal operation
    options.mode   = 'normal';
    P.mode = 'normal';
    powerDyn_SwingNormal = nonlinDASys(2,4,1,@SMIB_Swing,@SMIB_Swing_Const,options);
    options.tStart = 0.02; % Start time (tFinal of Fault)
    options.tFinal = 0.23; % Final time
    options.R0   = Rcont2{length(Rcont2)}{1}; % Initial state of last state from normal operation
    options.R0_y = Rcont_y2{length(Rcont_y2)}{1};
    [Rcont3,~,Rcont_y3] = reach(powerDyn_SwingNormal, options);
end


%% Random simulation results
runs = 30;
options.R0  = zonotope([options.x0,Bound_x]);
options.R_y = zonotope(options.y0.');

% M-matrix for ODE-15s solver
M = diag([ones(1,2) zeros(1,4)]);

for i=1:runs
    iChange = 1;
    while iChange <=3
        if iChange==1
            options.tStart = 0;
            options.tFinal = 0.01; 

            if i<=20
            options.x0=randPointExtreme(options.R0);
            else
            options.x0=randPoint(options.R0);
            end

        elseif iChange == 2
            options.tStart = 0.01; 
            options.tFinal = 0.02; 
            xFinal = x{i,1}(end,:)';
            options.x0 = xFinal(1:2);
            options.y0 = xFinal((2+1):end).';
        else
            options.tStart = 0.02;
            options.tFinal = 0.23;
            xFinal = x{i,2}(end,:)';
            options.x0 = xFinal(1:2);
            options.y0 = xFinal((2+1):end).';
        end

    options.M = odeset ('Mass',M,'RelTol',1e-8,'AbsTol',1e-10);
        [t{i,iChange},x{i,iChange}] = ode15s(@(t,x) SMIB_Model(t,x,iChange,I),[options.tStart, options.tFinal],[options.x0;options.y0.'],options.M); % Solve ODE
    iChange = iChange+1;
    end
end


%% Plotting
figure
hold on
projectedDimensions=[1 2];
for i=1:length(Rcont3)
    Zproj = project(Rcont3{i}{1},projectedDimensions);
    Zproj = reduce(Zproj,'girard',10);
    plotFilled(Zproj,[1 2],[.8 .8 .8],'EdgeColor','none');
end
for i=1:length(Rcont2)
    Zproj = project(Rcont2{i}{1},projectedDimensions);
    Zproj = reduce(Zproj,'girard',10);
    plotFilled(Zproj,[1 2],[.2 .2 .2],'EdgeColor','none');
end
for i=1:length(Rcont1)
    Zproj = project(Rcont1{i}{1},projectedDimensions);
    Zproj = reduce(Zproj,'girard',10);
    plotFilled(Zproj,[1 2],[.6 .6 .6],'EdgeColor','none');
end


for j=1:size(t,1)
    for iChange = 1:3
        plot(x{j,iChange}(:,projectedDimensions(1)),x{j,iChange}(:,projectedDimensions(2)),'k');
    end
end

%example completed
completed = 1;

function X0=init_network(X, P)
%% Variables
theta1 = X(1);
q1     = X(2);

% Power flow for SMIB transmission network
X0=[P.v1*P.v2/P.xs*sin(theta1-P.theta2)-P.p1;... 
    P.v1*P.v1/P.xs-P.v1*P.v2/P.xs*cos(theta1-P.theta2)-q1;...
];

function X0=init_generator(X,P,I)
%% Variables
delta1  = X(1);
omega1  = X(2);
E1      = X(3);
Pm1     = X(4);

% Here insert the equations which describe the system
X0=[  P.omegaS*(omega1);...
(1/P.M1)*(Pm1-I.P1 - (P.D1*omega1));...
I.P1 - (E1*I.V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-I.theta1));...
I.Q1  + (E1*I.V1*abs(P.xd1)*sin(angle(P.xd1)+delta1-I.theta1)) - (I.V1^2*abs(P.xd1)*sin(angle(P.xd1)));...

];


function out = SMIB_Model(t,x,iChange,I)
global P

xs = P.xs;

if iChange == 2
    xs = P.xs3;
end
%Ef=u(1);

%% Assign variables
%% Assign variables
delta1   = x(1);        % rotor angle (diff. variable)
omega1   = x(2);        % rotor angular speed (diff. variable)
P1       = x(3);        % active power at bus 1 (alg. variable)
Q1       = x(4);        % reactive power at bus 1 (alg. variable)
Theta1   = x(5);       % phase angle at bus 1 (alg. variable)
V1       = x(6);       % voltage at bus 1 (alg. variable)

%% The following steps are needed for realization of the controller
Pe1 = (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-Theta1));

%% Assign diff. and alg. variables
%  Differential variables
dx    = zeros(2,1);
dx(1) = P.omegaS*omega1;
dx(2) = (1/P.M1)*(I.Pm1-Pe1 - (P.D1*omega1));


%  Algebraic variables
dy    = zeros(4,1);
dy(1) = P1 - Pe1;
dy(2) = Q1 + (I.E1*V1*abs(P.xd1)*sin(angle(P.xd1)+delta1-Theta1)) - (V1^2*abs(P.xd1)*sin(angle(P.xd1)));
dy(3) = V1*P.v2/xs*sin(Theta1-P.theta2)-P1;
dy(4) = (V1*V1/xs)-V1*P.v2/xs*cos(Theta1-P.theta2)-Q1;

out = [dx;dy];
