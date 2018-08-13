function completed = example_reachLPV_SMIB()
% updated: 30-August-2017, AElG
% Implements reachability algorithm as in PES GE meeting 2017 paper:
% A. El-Guindy "LPV control with formal gurantees for transient stability
% of power systems"

% Define a global variable for all m.file 
global P

% SMIB system with 3 diff. variables and 8 alg. variables
dim_x = 3;
dim_y = 8;

reachflag = 1;

% Set options -------------------------------------------------------------
options.tStart   = 0;      % Start time
options.timeStep = 0.005;   % Time step size for reachable set computation

%% Assign fixed parameters ------------------------------------------------
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

%% Network parameters -----------------------------------------------------
P.xT        = 0.15;     % Reaktance of the transformers
P.xl1       = 0.5;      % Reaktance of the 1. line
P.xl2       = 0.93;     % Reaktance of the 2. line

%% Resulting overall network reactance for normal and fault case ----------
P.xs        = P.xT + 1/(1/P.xl1 + 1/P.xl2); % Pre-fault
P.xs2       = P.xT; % Fault
P.xs3       = P.xT + P.xl1; % Post-fault

%% Power Flow: Initial values ---------------------------------------------
% Bus-1 is PV-bus
P.v1        = 1;           % voltage at bus-1
P.p1        = 0.9;         % active power at bus-1

% Bus-2 is slack bus
P.v2        = 0.90081;     % Voltage at Infinite Bus := const
P.theta2    = 0;           % Phase at Infinite Bus := const

% Set options for nonlinear solver to obtain initial variables of remaining
% variables of the system
options_fsolve = optimoptions('fsolve','Display','off');
X0netz = fsolve(@(X)init_sim_netzwerk(X, P), ones(4,1),options_fsolve);
% Zuweisung der Variablen
P.theta1  = X0netz(1);       % phase at bus-1
P.q1      = X0netz(2);       % reactive power at bus-1

%% Generator Initinal values ----------------------------------------------
X0gen=fsolve(@(X)init_generator(X, P), ones(10,1),options_fsolve);

% diff. variables
delta0  = X0gen(1);        % Zustand - state variable
omega0  = X0gen(2);        % Zustand
Eq_0    = X0gen(3);        % Zustand

x0 = [ delta0 omega0 Eq_0];

% alg. variables
ed0     = X0gen(4);        % alg. Variable
eq0     = X0gen(5);        % alg. Variable
Ef      = X0gen(6);        % Constant field voltage
id0     = X0gen(7);        % alg. Variable
iq0     = X0gen(8);        % alg. Variable
P.Tm    = X0gen(9);        % Constant mechanical torque
P.Te    = X0gen(10);       % alg. Variable

y0 = [ ed0 eq0 id0 iq0 P.p1 P.q1 P.theta1 P.v1];

%% Set value of diff. variables
u0 = [1.2242; 0 ;1.0241];

%% Assign parameters for reachability analysis ----------------------------


options.x0 = x0;
options.y0 = y0; 
options.y0guess = y0.';
Bound_x = [0.005  0 0; 0 0.5e-4 0; 0 0 0.005];

options.R0=zonotope([options.x0.',Bound_x]); % x0
options.R_y=zonotope([options.y0.',zeros(dim_y,dim_y)]); % y0

%% Initial set for the input variables ------------------------------------
options.uTrans = u0;
Bound_u(3,3)=0;                           % Uncertainty
options.U=zonotope([zeros(3,1),Bound_u]); % Initial input state
options.originContained = 0;              % No origin

options.zonotopeOrder = 300; % Zonotope order
options.polytopeOrder = 10;  % Polytope order

options.taylorTerms = 4;

options.tensorOrder          = 2;
options.advancedLinErrorComp = 1;
options.intermediateOrder    = 3;

options.reductionTechnique = 'girard';

options.errorOrder = 3;
options.maxError   = 0;
options.oldError_x = 0.01*ones(dim_x,1);
options.maxError_x = options.oldError_x;
options.oldError_y = 0.005*ones(dim_y,1);
options.maxError_y = options.oldError_y;

for n=1:(options.taylorTerms+1)
    %time step
    r = options.timeStep;
    %compute initial state factor
    options.factor(n)= r^(n)/factorial(n);    
end

options.reductionInterval = 1e5;
%--------------------------------------------------------------------------

% Specify continuous dynamics----------------------------------------------
% I have done this step already !!
%powerDyn = nonlinDASys(3,8,3,@SMIB_States_Paper,@SMIB_Const_Paper,options);
%save powerDyn powerDyn
%load powerDynFault
%load powerDyn

%--------------------------------------------------------------------------

% Compute reachable set
if reachflag
%% normal operation
    P.mode       = 'normal';
    powerDyn = nonlinDASys(3,8,3,@SMIB_States_Paper,@SMIB_Const_Paper,options);
    options.tStart = 0;   % Start time
    options.tFinal = 0.1; % Final time
    [Rcont1,~,Rcont_y1] = reach(powerDyn, options);
    
    %% Fault
    
    P.mode       = 'fault';
    powerDynFault = nonlinDASys(3,8,3,@SMIB_States_Paper,@SMIB_Const_Paper,options);
    options.tStart = 0.1; % Start time (tFinal of normal operation)
    options.tFinal = 0.3; % Final time
    options.R0   = Rcont1{length(Rcont1)}{1}; % Initial state of last state from normal operation
    options.R0_y = Rcont_y1{length(Rcont_y1)}{1};
    [Rcont2,~,Rcont_y2] = reach(powerDynFault, options);

    %% normal operation
    P.mode = 'normal';
    powerDyn = nonlinDASys(3,8,3,@SMIB_States_Paper,@SMIB_Const_Paper,options);
    options.tStart = 0.3; % Start time (tFinal of Fault)
    options.tFinal = 1; % Final time
    options.R0   = Rcont2{length(Rcont2)}{1}; % Initial state of last state from normal operation
    options.R0_y = Rcont_y2{length(Rcont_y2)}{1};
    [Rcont3,~,Rcont_y3] = reach(powerDyn, options);
    
end

for j=1:length(Rcont1)
    Xhull1{j}        = interval(Rcont1{j}{1});
    delta_para1{j}   = interval(infimum(Xhull1{j}(1)),supremum(Xhull1{j}(1)));
    Yhull1{j}        = interval(Rcont_y1{j}{1});
    id1{j}           = interval(infimum(Yhull1{j}(3)),supremum(Yhull1{j}(3)));
    iq1{j}           = interval(infimum(Yhull1{j}(4)),supremum(Yhull1{j}(4)));
    
    theata11{j} = (P.Tm-(P.xq-P.x_d)*id1{j}*iq1{j})/(2*P.H*delta_para1{j} );
    theata21{j} = iq1{j};
    theata31{j} = id1{j}/delta_para1{j};
end

for j=1:length(Rcont2)
    Xhull2{j}        = interval(Rcont2{j}{1});
    delta_para2{j}   = interval(infimum(Xhull2{j}(1)),supremum(Xhull2{j}(1)));
    Yhull2{j}        = interval(Rcont_y2{j}{1});
    id2{j}           = interval(infimum(Yhull2{j}(3)),supremum(Yhull2{j}(3)));
    iq2{j}           = interval(infimum(Yhull2{j}(4)),supremum(Yhull2{j}(4)));
    
    
    theata12{j} = (P.Tm-(P.xq-P.x_d)*id2{j}*iq2{j})/(2*P.H*delta_para2{j} );
    theata22{j} = iq2{j};
    theata32{j} = id2{j}/delta_para2{j};
end

for j=1:length(Rcont3)
    Xhull3{j}        = interval(Rcont3{j}{1});
    delta_para3{j}   = interval(infimum(Xhull3{j}(1)),supremum(Xhull3{j}(1)));
    Yhull3{j}        = interval(Rcont_y3{j}{1});
    id3{j}           = interval(infimum(Yhull3{j}(3)),supremum(Yhull3{j}(3)));
    iq3{j}           = interval(infimum(Yhull3{j}(4)),supremum(Yhull3{j}(4)));
    
    theata13{j} = (P.Tm-(P.xq-P.x_d)*id3{j}*iq3{j})/(2*P.H*delta_para3{j} );
    theata23{j} = iq3{j};
    theata33{j} = id3{j}/delta_para3{j};
end

theata1 = [theata11,theata12,theata13];
theata2 = [theata21,theata22,theata23];
theata3 = [theata31,theata32,theata33];

figure
subplot(3,1,1)
hold on
for i=1:length(theata1)
    IH1 = interval([infimum(theata1{i}); infimum(theata1{i})], [supremum(theata1{i});supremum(theata1{i})]);
    t1 = (i-1)*options.timeStep;
    t2 = i*options.timeStep;
    IH1 = interval([t1;infimum(IH1(1))], [t2;supremum(IH1(1))]);
    plotFilled(IH1,[1 2],[.2 .2 .2],'EdgeColor','none');
end

subplot(3,1,2)
hold on
for i=1:length(theata2)
    IH1 = interval([infimum(theata2{i}); infimum(theata2{i})], [supremum(theata2{i});supremum(theata2{i})]);
    t1 = (i-1)*options.timeStep;
    t2 = i*options.timeStep;
    IH1 = interval([t1;infimum(IH1(1))], [t2;supremum(IH1(1))]);
    plotFilled(IH1,[1 2],[.2 .2 .2],'EdgeColor','none');
end

subplot(3,1,3)
hold on
for i=1:length(theata3)
    IH1 = interval([infimum(theata3{i}); infimum(theata3{i})], [supremum(theata3{i});supremum(theata3{i})]);
    t1 = (i-1)*options.timeStep;
    t2 = i*options.timeStep;
    IH1 = interval([t1;infimum(IH1(1))], [t2;supremum(IH1(1))]);
    plotFilled(IH1,[1 2],[.2 .2 .2],'EdgeColor','none');
end

%% Obtain interval hull of all reachable sets Eq. (13) in paper
% Uncomment after implementation of enclose
% for j=1:length(Rcont1)
%     XInt1{j} = interval(Rcont1{j}{1});
%     YInt1{j} = interval(Rcont_y1{j}{1});
% end
% 
% for j=1:length(Rcont2)
%     XInt2{j} = interval(Rcont2{j}{1});
%     YInt2{j} = interval(Rcont_y2{j}{1});
% end
% 
% for j=1:length(Rcont3)
%     XInt3{j} = interval(Rcont3{j}{1});
%     YInt3{j} = interval(Rcont_y3{j}{1});
% end
% 
% XInt=[XInt1,XInt2,XInt3];
% YInt=[YInt1,YInt2,YInt3];

%% Enclose all hulls into one interval Eq. (14)

% X_states = enclose(XInt{1},XInt{2});
% Y_states = enclose(YInt{1},YInt{2});
% 
% for i=3:length(YInt)
%     X_states = enclose(X_states,YInt{i});
%     Y_states = enclose(Y_states,YInt{i});
% end
% 
% % Obtain upper and lower bound of each variable
% 
% delta_para = interval(infimum(X_states(1)),supremum(X_states(1)));
% id_para    = interval(infimum(Y_states(3)),supremum(Y_states(3)));
% id_para    = interval(infimum(Y_states(4)),supremum(Y_states(4)));
% 
% % Obtain bounds of the scheduling parameters \theta
% theata1_para = (P.Tm-(P.xq-P.x_d)*id_para*iq_para)/(2*P.H*delta_para );
% theata2_para = iq_para;
% theata3_para = id_para/delta_para;

%% Random simulation results
runs = 30;
options.R0=zonotope([options.x0.',Bound_x]);
options.R_y=zonotope(options.y0.');
options.R0 = reduce(options.R0,'girard',10);

% M-matrix for ODE-15s solver
M = diag([ones(1,3) zeros(1,8)]);

for i=1:runs
    iChange = 1;
    while iChange <=3
        if iChange==1
            options.tStart = 0;
            options.tFinal = 0.1; 

            if i<=20
            options.x0=randPointExtreme(options.R0);
            else
            options.x0=randPoint(options.R0);
            end

        elseif iChange == 2
            options.tStart = 0.1; 
            options.tFinal = 0.3; 
            xFinal = x{i,1}(end,:)';
            options.x0 = xFinal(1:3);
            options.y0 = xFinal((3+1):end).';
        else
            options.tStart = 0.3;
            options.tFinal = 1;
            xFinal = x{i,2}(end,:)';
            options.x0 = xFinal(1:3);
            options.y0 = xFinal((3+1):end).';
        end

    
   % Input
    if i<=20
        options.u=randPointExtreme(options.U)+options.uTrans;
    else
        options.u=randPoint(options.U)+options.uTrans;
    end

    options.M = odeset ('Mass',M,'RelTol',1e-8,'AbsTol',1e-10);
        [t{i,iChange},x{i,iChange}] = ode15s(@(t,x) SMIB_Model(t,x,iChange,P),[options.tStart, options.tFinal],[options.x0;options.y0.'],options.M); % Solve ODE
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
        plot(x{j,iChange}(:,projectedDimensions(1)),x{j,iChange}(:,projectedDimensions(2)),'Color',iChange*[0.01 0.1 0.3]);
    end
end
plotFilled(options.R0,[1 2],[1 1 1],'EdgeColor',[.2 .2 .2]); 

%example completed
completed = 1;

function X0=init_sim_netzwerk(X, P)
%% Initialisation of the grid, POWER FLOW
theta1 = X(1);
q1     = X(2);
p2     = X(3);
q2     = X(4);

X0=[P.v1*P.v2/P.xs*sin(theta1-P.theta2)-P.p1;...       %Netzwerk
    P.v1*P.v2/P.xs*sin(P.theta2-theta1)-p2;...
    P.v1*P.v1/P.xs-P.v1*P.v2/P.xs*cos(theta1-P.theta2)-q1;...
    P.v2*P.v2/P.xs-P.v1*P.v2/P.xs*cos(P.theta2-theta1)-q2;...
];

function X0=init_generator(X,P)
%% Assignment of the variables
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

function out = SMIB_Model(t,x,iChange,P)
%global P

xs = P.xs;

if iChange == 2
    xs = P.xs3;
end
%Ef=u(1);

%% Assign variables
delta   = x(1);        %Zustand
omega   = x(2);        %Zustand
Eq_     = x(3);        %Zustand 
ed      = x(4);        %alg. Variable
eq      = x(5);        %alg. Variable
id      = x(6);        %alg. Variable
iq      = x(7);        %alg. Variable
P1      = x(8);        %alg. Variable
Q1      = x(9);        %alg. Variable
Theta1  = x(10);       %alg. Variable
V1      = x(11);       %alg. Variable
%ureal   = x(12);       %alg. Variable

%% Controller gain
x1 = (P.Tm-(P.xq-P.x_d)*id*iq)/(2*P.H*delta); 
x2 = iq; 
x3 = id/delta;
  
%% Parameter Depdandant Controller ----------------------------
% See Eq. (9)
% Obtained using the function closed-form expressions
alpha = [ -(25*x2 - 391/40)*((500*x1)/17 - 121/68)*((500*x3)/133 - 881/266);
           (25*x2 - 391/40)*((500*x1)/17 - 121/68)*((500*x3)/133 - 615/266);
           (25*x2 - 351/40)*((500*x1)/17 - 121/68)*((500*x3)/133 - 881/266);
          -(25*x2 - 351/40)*((500*x1)/17 - 121/68)*((500*x3)/133 - 615/266);
           (25*x2 - 391/40)*((500*x1)/17 - 53/68)*((500*x3)/133 - 881/266);
          -(25*x2 - 391/40)*((500*x1)/17 - 53/68)*((500*x3)/133 - 615/266);
          -(25*x2 - 351/40)*((500*x1)/17 - 53/68)*((500*x3)/133 - 881/266);
           (25*x2 - 351/40)*((500*x1)/17 - 53/68)*((500*x3)/133 - 615/266) ];
       
K = 1.0e+04 * [0.011422432372455   1.509509521267021  -0.011003653524966;
               0.011462598370783   1.509509521189578  -0.011003653524654;
               0.011190584090320   1.447666229394653  -0.011358761159364;
               0.011230716616862   1.447664877138313  -0.011358755674191;
               0.018119719099198   1.606155453002849  -0.011331758629217;
               0.018159885100296   1.606155452990147  -0.011331758629685;
               0.017383029200531   1.567306261121143  -0.011736512340165;
               0.017423180141780   1.567308539700949  -0.011736540258299];

%% Control error
deltarl = 1.2242;
omegarl = 0;
Eq_rl   = 1.0241;

Kp = alpha.'*K;
vf = Kp*[delta-deltarl; omega-omegarl ;Eq_-Eq_rl]+2.4207;


%% State Varaibles
dx(1) = P.omegaS*omega;
dx(2) = ((P.Tm-(P.xq-P.x_d)*id*iq)-iq*Eq_-P.D*omega)/(2*P.H);
dx(3) = (vf-(P.xd-P.x_d)*id-Eq_)/P.T_d0;

%% Alg. Variables;
dx(4) = -ed+V1*sin(delta-Theta1);
dx(5) = -eq+V1*cos(delta-Theta1);
dx(6) = -ed+P.xq*iq-P.ra*id;
dx(7) = -eq+(Eq_)-P.ra*iq-P.x_d*id;
dx(8) = -P1+ed*id+eq*iq;
dx(9) = -Q1+eq*id-ed*iq;
dx(10) = V1*P.v2/xs*sin(Theta1-P.theta2)-P1;
dx(11) = V1*V1/xs-V1*P.v2/xs*cos(Theta1-P.theta2)-Q1;

out = [dx].';
