function completed = example_composite_9bus()
% updated: 30-August-2017, AElG
% As realized in the ACC paper
% Compositional transient stability anaylsis of Power Systems

clear
clc

global P I
[P,~,I]   = GetStarted();

options = SetOptions(2);


% 9-bus system with 3 diff. variables and 8 alg. variables
y_gen = [ I.delta1 I.omega1 ....
          I.delta2 I.omega2 ...
          I.delta3 I.omega3 ];
        
y_grid = [  I.theta2 I.theta3 I.theta4 I.theta5 I.theta6 I.theta7 I.theta8 I.theta9...
            I.V1 I.V2 I.V3 I.v4 I.v5 I.v6 I.v7 I.v8 I.v9 ...
            I.Pe1 I.Pe2 I.Pe3 I.p4 I.p5 I.p6 I.p7 I.p8 I.p9...
            I.Q1 I.Q2 I.Q3 I.q4 I.q5 I.q6 I.q7 I.q8 I.q9];
            
  
% %% Assign parameters for reachability analysis ----------------------------
x0{1} = y_gen(1,1:2);
x0{2} = y_gen(1,3:4);
x0{3} = y_gen(1,5:6);

g_R0  = 1e-3*eye(6,6);

i = 1;
%--------------------------------------------------------------------------
%% Set classes
powerDyn_Gen1 = nonlinearSys(2,1,@Gen1_LPV,options);
powerDyn_Gen2 = nonlinearSys(2,2,@Gen2_LPV,options);
powerDyn_Gen3 = nonlinearSys(2,2,@Gen3_LPV,options);

t                = 0;

Rnext1 = zonotope([x0{1}.',g_R0(1:2,:)]);
Z01    = Rnext1;
Rnext2 = zonotope([x0{2}.',g_R0(3:4,:)]);
Z02    = Rnext2;
Rnext3 = zonotope([x0{3}.',g_R0(5:6,:)]); % x0
Z03    = Rnext3;

options1 = options;
options2 = options;
options3 = options;
options_fsolve = optimoptions('fsolve','Display','off');

tic
if options.reachflag
%% normal operation
    while t < options.tFinal
        
        % Check fault conditions
        if t > 0.01 && t < 0.05
            options.mode = 'fault';
            G.xs = P.xs1;
        else
            options.mode = 'normal';
            G.xs = P.xs;
        end
        
        % Increment time and set counter
        y_gen1   = center(Rnext1);
        G.delta1 = y_gen1(1);

        y_gen2   = center(Rnext2);
        G.delta2 = y_gen2(1);

        y_gen3   = center(Rnext3);
        G.delta3 = y_gen3(1);

        init_grid = @(X)PowerGrid(X,G,P,I);
        X0        = fsolve(init_grid, [zeros(8,1);ones(27,1)],options_fsolve );
        
        Rnext    = cartesianProduct(Rnext1,Rnext2);
        Rnext    = cartesianProduct(Rnext,Rnext3);
        
            delta_x       = center(Rnext);
            Const_gx_gy   = PowerGrid(X0,G,P,I);
            if strcmp(options.mode,'normal')
                [w_bar,J_bar] = jacobian_PowerGridNormal(X0,delta_x);
            elseif strcmp(options.mode,'fault')
                [w_bar,J_bar] = jacobian_PowerGridFault(X0,delta_x);
            end
            U_const = -Const_gx_gy + J_bar*delta_x + w_bar*X0;
            R_y{i} = pinv(w_bar)*(U_const + (-1*J_bar*Rnext));
            R_y{i} = deleteZeros(R_y{i});
        
        % Gen1
        options1.R0        = Rnext1; % x0
        options1.uTrans    = X0(9,1);
        options1.U         = zonotope([0 0.0045]);
        
        [Rcont1,Rnext1,options1] = contReach(powerDyn_Gen1,options1);

        % Gen2
        options2.R0        = Rnext2; % x0
        options2.uTrans    = [X0(1,1);X0(10,1)];
        options2.U         = zonotope([zeros(2,1) [0.00035 0;0 0.0025]]);

        [Rcont2,Rnext2,options2] = contReach(powerDyn_Gen2,options2);

        % Gen3
        options3.R0        = Rnext3; % x0
        options3.uTrans    = [X0(2,1);X0(11,1)];
        options3.U         = zonotope([zeros(2,1) [0.00035 0;0 0.0025]]);
        
        [Rcont3,Rnext3,options3] = contReach(powerDyn_Gen3,options3);
        
        Rcont{i}    = cartesianProduct(Rcont1,Rcont2);
        Rcont{i}    = cartesianProduct(Rcont{i},Rcont3);
        
        subplot(1,3,1)
        hold on
        plotFilled(Rcont{i},[1 2],[.8 .8 .8],'EdgeColor','none'); 
        subplot(1,3,2)
        hold on
        plotFilled(Rcont{i},[3 4],[.8 .8 .8],'EdgeColor','none'); 
        subplot(1,3,3)
        hold on
        plotFilled(Rcont{i},[5 6],[.8 .8 .8],'EdgeColor','none'); 
        
         
        t = t + options.timeStep;
        
         
       
        i = i+1;

    end
    
end
toc

runs = 10;
for i=1:runs
    %set initial state, input
    if i<=5
        y_gen1=randPointExtreme(Z01); % At the borders
        y_gen2=randPointExtreme(Z02); % At the borders
        y_gen3=randPointExtreme(Z03); % At the borders
    else
        y_gen1=randPoint(Z01); % Inside the Zonotope
        y_gen2=randPoint(Z02); % Inside the Zonotope
        y_gen3=randPoint(Z03); % Inside the Zonotope
    end

    
y_grid = [  I.theta2 I.theta3 I.theta4 I.theta5 I.theta6 I.theta7 I.theta8 I.theta9...
            I.V1 I.V2 I.V3 I.v4 I.v5 I.v6 I.v7 I.v8 I.v9 ...
            I.Pe1 I.Pe2 I.Pe3 I.p4 I.p5 I.p6 I.p7 I.p8 I.p9...
            I.Q1 I.Q2 I.Q3 I.q4 I.q5 I.q6 I.q7 I.q8 I.q9];
            
y_gen = [y_gen1;y_gen2;y_gen3].'  ;  
y0 = [ y_gen y_grid].';

solveDAE(y_gen,y0,P,I)
end

subplot(1,3,1)
hold on
plotFilled(Z01,[1 2],[1 1 1],'EdgeColor',[.2 .2 .2]); 
subplot(1,3,2)
hold on
plotFilled(Z02,[1 2],[1 1 1],'EdgeColor',[.2 .2 .2]); 
subplot(1,3,3)
hold on
plotFilled(Z03,[1 2],[1 1 1],'EdgeColor',[.2 .2 .2]); 

%example completed
completed = 1;


function X0_G = PowerGrid(x,G,P,I)
xs      = G.xs;
delta1  = G.delta1;
delta2  = G.delta2;
delta3  = G.delta3;

%% Alg. State variable of the Power Flow
% Bus phase angle
theata1 = I.theta1;
theata2 = x(1);
theata3 = x(2);
theata4 = x(3);
theata5 = x(4);
theata6 = x(5);
theata7 = x(6);
theata8 = x(7);
theata9 = x(8);

% Bus voltage levels
V1 = x(9);
V2 = x(10);
V3 = x(11);
v4 = x(12);
v5 = x(13);
v6 = x(14);
v7 = x(15);
v8 = x(16);
v9 = x(17);

p1 = x(18);
p2 = x(19);
p3 = x(20);
p4 = x(21);
p5 = x(22); 
p6 = x(23);
p7 = x(24);
p8 = x(25);
p9 = x(26);

q1 = x(27);
q2 = x(28);
q3 = x(29);
q4 = x(30);
q5 = x(31);
q6 = x(32);
q7 = x(33);
q8 = x(34);
q9 = x(35);

Pw    = [ p1 p2 p3 p4 p5 p6 p7 p8 p9]';
Q     = [ q1 q2 q3 q4 q5 q6 q7 q8 q9]';            
V     = [V1 V2 V3 v4 v5 v6 v7 v8 v9].';
Theta = [theata1 theata2 theata3 theata4 theata5 theata6 theata7 theata8 theata9];

n        = length (V);
Theta_hk = Theta.'*ones(1,n)-ones(n,1)*Theta ;

X0_G = [ p1 - (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-theata1))
         p2 - (I.E2*V2*abs(P.xd2)*cos(angle(P.xd2)+delta2-theata2));
         p3 - (I.E3*V3*abs(P.xd3)*cos(angle(P.xd3)+delta3-theata3));
p4 - I.p4 * v4^2 / I.v4^2
p5 - I.p5 * v5^2 / I.v5^2
p6 - I.p6 * v6^2 / I.v6^2
p7 - I.p7 * v7^2 / I.v7^2
p8 - I.p8 * v8^2 / I.v8^2
p9 - I.p9 * v9^2 / I.v9^2

%dx(10)  = q1 - (I.E1*V1*abs(P.xd1)*sin(angle(P.xd1)+delta1-theata1)) + (V1^2*abs(P.xd1)*sin(angle(P.xd1)));
q2 + (I.E2*V2*abs(P.xd2)*sin(angle(P.xd2)+delta2-theata2)) - (V2^2*abs(P.xd2)*sin(angle(P.xd2)))
q3 + (I.E3*V3*abs(P.xd3)*sin(angle(P.xd3)+delta3-theata3)) - (V3^2*abs(P.xd3)*sin(angle(P.xd3)))
q4 - I.q4 * v4^2 / I.v4^2
q5 - I.q5 * v5^2 / I.v5^2
q6 - I.q6 * v6^2 / I.v6^2
q7 - I.q7 * v7^2 / I.v7^2
q8 - I.q8 * v8^2 / I.v8^2
q9 - I.q9 * v9^2 / I.v9^2

V.*( (real(xs).*cos(Theta_hk) + imag(xs).*sin(Theta_hk))* V ) - Pw 
V.*( (real(xs).*sin(Theta_hk) - imag(xs).*cos(Theta_hk))* V ) - Q 
];
  
function [Rcont, Rnext,options] = contReach(sys,options)
% initialize reachable set computations
[Rnext_x,options] = initReach(sys, options.R0, options);

% save reachable set in cell structure
Rnext     = reduce(Rnext_x.tp{1}.set, 'girard', 5); 
Rcont     = reduce(Rnext_x.ti{1}, 'girard', 5); 
  
function solveDAE(y_gen,y0,P,I) 
t1 = 0;
t2 = 0.01;
t3 = 0.05;
t4 = 0.3;
tspan1 = [t1 t2];
tspan2 = [t2 t3];
tspan3 = [t3 t4];

%% Solve using ODE-15s, for the transient stablity sequence
%  Define M-matrix for diff. and alg. varaibles
DGL_anz = length(y_gen);    % L??nge der DGL f??r Mass- Matrix
ALG_anz = length(y0) - DGL_anz;                         % L??nge der ALG f??r Mass- Matrix

M = diag([ones(1,DGL_anz) zeros(1,ALG_anz)]);
options_fsolve = odeset ('Mass',M,'RelTol',1e-05,'MaxStep',0.0001);

xs = P.xs;
f1r = @(t,x)LPV_3Gen(t,x,xs,I,P);
xs = P.xs1;
f2r = @(t,x)LPV_3Gen(t,x,xs,I,P);
xs = P.xs;
f3r = @(t,x)LPV_3Gen(t,x,xs,I,P);

[~,y1r] = ode15s(f1r,tspan1,y0,options_fsolve);
[~,y2r] = ode15s(f2r,tspan2,y1r(end,:),options_fsolve);
[~,y3r] = ode15s(f3r,tspan3,y2r(end,:),options_fsolve);

%Tr=[t1r;t2r;t3r];

Yr=[y1r;y2r;y3r];

%Ef=Yr(:,[61, 62,63]);


figure(1)
hold on
subplot(1,3,1)
 plot(Yr(:,1), Yr(:,2),'k')
 subplot(1,3,2)
 plot(Yr(:,3), Yr(:,4),'k')
 subplot(1,3,3)
 plot(Yr(:,5),Yr(:,6),'k')
 

function [options] = SetOptions(dim_x)
%global P

options.Type = 'Power System';
options.reachflag = 1;

options.tStart   = 0;   % Start time
options.tFinal   = 0.13;   % Final time 
options.timeStep = 0.0005;

% Set options -------------------------------------------------------------
% Initial set for the input variables ------------------------------------
options.originContained = 0;              % No origin

options.zonotopeOrder = 300; % Zonotope order
options.polytopeOrder = 10;  % Polytope order

options.taylorTerms = 15;

options.tensorOrder          = 2;
options.advancedLinErrorComp = 1;
options.intermediateOrder    = 3;

options.reductionTechnique = 'girard';

options.errorOrder = 3;
options.oldError   = 1e-3*ones(dim_x,1);
options.maxError   = 100*options.oldError;

options.reductionInterval = 1e5;

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %time step
    r = options.timeStep;
    %compute initial state factor
    options.factor(i)= r^(i)/factorial(i);    
end

function [P,G,I] = GetStarted()
%% Assign fixed parameters of the synchronous generator
% 1. Generator
P.M1         = 1/(15*pi);   	% Tr???gheitskoeffizient
P.D1         = 0.8;           % D???mpfungskoeffizient
P.xd1       = 1/(1i*0.2);  	% d-Achsen transiente Reaktanz     
P.omegaS     = 2*pi*50;

% 2. Generator
P.M2         = 1/(15*pi);   	 % Tr???gheitskoeffizient
P.D2         = 0.8;        % D???mpfungskoeffizient
P.xd2       = 1/(1i*0.2);  	% d-Achsen transiente Reaktanz

% 3. Generator
P.M3         = 1/(15*pi);   	% Tr???gheitskoeffizient
P.D3         = 0.8;        % D???mpfungskoeffizient
P.xd3       = 1/(1i*0.2);  	% d-Achsen transiente Reaktanz

%% Network parameters
% Transformers
G.xT1        =  1i*0.0576;                                        % Reaktanz des 1. Trafo
G.xT2        =  1i*0.0625;                                        % Reaktanz des 2. Trafo
G.xT3        =  1i*0.0586;                                        % Reaktanz des 3. Trafo

% Bus-impedances
G.xl45 = 0.01+1i*0.085;     % Impedanz Bus 4-5
G.xl46 = 0.017+1i*0.092;    % Impedanz Bus 4-6
G.xl57 = 0.032 +1i*0.161;   % Impedanz Bus 5-7
G.xl69 = 0.039+1i*0.170;    % Impedanz Bus 6-9
G.xl78 = 0.0085+1i*0.072;   % Impedanz Bus 7-8
G.xl89 = 0.0119+1i*0.1008;  % Impedanz Bus 8-9

% bus susceptances
G.b45 = 0.088*1i;               % Suszeptanz der einzelnen Leitungen
G.b46 = 0.079*1i;
G.b57 = 0.153*1i;
G.b69 = 0.179*1i;
G.b78 = 0.0745*1i;
G.b89 = 0.1045*1i;

% Definition of the ADMITTANCE MATRIX --- Aufstellen der Admittanzmatrix f???r den Power Flow
% Original-System as describes by Anderson and FOuad --- Originalsystem
G.Y         = [1/G.xT1             0               0               -1/G.xT1               0             0               0              0               0;
                0               1/G.xT2              0                   0                  0             0              -1/G.xT2        0               0;
                0                   0             1/G.xT3                0                  0             0               0              0          -1/G.xT3;
               -1/G.xT1               0               0  (1/G.xT1+1/G.xl45+1/G.xl46+G.b45+G.b46) -1/G.xl45  -1/G.xl46       0              0               0;
                0                   0               0 -1/G.xl45 (1/G.xl45+1/G.xl57+G.b45+G.b57)          0           -1/G.xl57           0              0;
                0                   0               0 -1/G.xl46         0       (1/G.xl46+1/G.xl69+G.b46+G.b69)            0            0            -1/G.xl69;
                0           -1/G.xT2        0       0 -1/G.xl57         0                   (1/G.xT2+1/G.xl57+1/G.xl78+G.b57+G.b78)          -1/G.xl78  0;
                0           0            0      0           0       0           -1/G.xl78               (1/G.xl78+1/G.xl89+G.b78+G.b89)         -1/G.xl89;
                0           0        -1/G.xT3       0       0       -1/G.xl69       0                       -1/G.xl89  (1/G.xT3+1/G.xl69+1/G.xl89+G.b89+G.b69);
                ];

%% Admittance-matrizes for different faults
% Admittance change of line from Bus 5 to Bus 7: half of the admittances
G.Yf         = [ 1/G.xT1              0               0               -1/G.xT1               0             0               0              0               0;
                0               1/G.xT2              0                   0                  0             0              -1/G.xT2        0               0;
                0                   0             1/G.xT3                0                  0             0               0              0          -1/G.xT3;
             -1/G.xT1               0               0  (1/G.xT1+1/G.xl45+1/G.xl46+G.b45+G.b46) -1/G.xl45  -1/G.xl46       0              0               0;
                0                   0               0 -1/G.xl45 (1/G.xl45+2/G.xl57+G.b45+2*G.b57)          0           -2/G.xl57           0              0;
                0                   0               0 -1/G.xl46         0       (1/G.xl46+1/G.xl69+G.b46+G.b69)            0            0            -1/G.xl69;
                0           -1/G.xT2        0       0 -2/G.xl57         0                   (1/G.xT2+2/G.xl57+1/G.xl78+2*G.b57+G.b78)          -1/G.xl78  0;
                0           0            0      0           0       0           -1/G.xl78               (1/G.xl78+1/G.xl89+G.b78+G.b89)         -1/G.xl89;
                0           0        -1/G.xT3       0       0       -1/G.xl69       0                       -1/G.xl89  (1/G.xT3+1/G.xl69+1/G.xl89+G.b89+G.b69);
                ];


xs = G.Yf;        % Originale Admittanzmatrix 
G.Ystart=xs;
xs1 = G.Y;        % Admittanzmatrix mit Fehler G.Yf=nur Halbierung der Leitung, G.Yf1=Schlimm

P.xs  = xs;  % Normal operation
P.xs1 = xs1; % Fault


%% Power Flow: Startvorgaben
G.v1         = 1.04;         % Spannung am Slack Bus (Referenzknoten 1)
G.theta1    = 0;             % Winkel am Slack Bus
G.p2        = 1.63;          % Wirkleistung an Gen. 2
G.v2        = 1.025;         % Spannung an Gen. 2
G.p3        = 0.85;          % Wirkleistung an Gen. 3
G.v3        = 1.025;         % Spannung an Gen. 3
G.p4        = 0;
G.p5        = -1.25;
G.p6        = -0.5;
G.p7        = 0;
G.p8        = -1;
G.p9        = 0;
G.q4        = 0;
G.q5        = -0.5;
G.q6        = -0.2842;
G.q7        = 0;
G.q8        = -0.35;
G.q9        = 0;

%Solution of the Power Flows Problem
options = optimoptions('fsolve','Display','off');

init_grid = @(X)A_init_grid_9k(X,G);
X0 = fsolve(init_grid, ones(18,1),options);


% Initial Values of the grid after Power Flow
I.Pe1 = X0(1);
I.Pe2 = G.p2;
I.Pe3 = G.p3;
I.p4 = G.p4;
I.p5 = G.p5;
I.p6 = G.p6;
I.p7 = G.p7; 
I.p8 = G.p8;
I.p9 = G.p9;

I.Q1 = X0(2);
I.Q2 = X0(3);
I.Q3 = X0(4);
I.q4 = G.q4;
I.q5 = G.q5;
I.q6 = G.q6;
I.q7 = G.q7; 
I.q8 = G.q8;
I.q9 = G.q9;

I.theta1 = G.theta1;
I.theta2 = X0(5);
I.theta3 = X0(6);
I.theta4 = X0(7);
I.theta5 = X0(8);
I.theta6 = X0(9);
I.theta7 = X0(10); 
I.theta8 = X0(11);
I.theta9 = X0(12);

I.V1 = G.v1;
I.V2 = G.v2;
I.V3 = G.v3;
I.v4 = X0(13);
I.v5 = X0(14);
I.v6 = X0(15);
I.v7 = X0(16); 
I.v8 = X0(17);
I.v9 = X0(18);

alg = @(X)A_init_gen(X,P,I);
X0_gen = fsolve(alg, ones(12,1),options);

% Generator 1
I.delta1 = X0_gen(1);
I.omega1 = 0;
I.Pm1    = X0_gen(3);
I.E1     = X0_gen(4);

% Generator 2
I.delta2 = X0_gen(5);
I.omega2 = 0;
I.Pm2    = X0_gen(7);
I.E2     = X0_gen(8);

% Generator 3
I.delta3 = X0_gen(9);
I.omega3 = 0;
I.Pm3    = X0_gen(11);
I.E3     = X0_gen(12);

function X0_G = A_init_grid_9k (X,G)

% Parameter f???r das Netz
p1 = X(1);
p2 = G.p2;
p3 = G.p3;
p4 = G.p4;
p5 = G.p5;
p6 = G.p6;
p7 = G.p7;  % X(13)
p8 = G.p8;
p9 = G.p9;

q1 = X(2);
q2 = X(3);
q3 = X(4);
q4 = G.q4;
q5 = G.q5;
q6 = G.q6;
q7 = G.q7;  % ;X(14)
q8 = G.q8;
q9 = G.q9;

theta1 = G.theta1;
theta2 = X(5);
theta3 = X(6);
theta4 = X(7);
theta5 = X(8);
theta6 = X(9);
theta7 = X(10); %0.21431 Wert aus PSAT (8. Knoten entnommen)
theta8 = X(11);
theta9 = X(12);

v1 = G.v1;
v2 = G.v2;
v3 = G.v3;
v4 = X(13);
v5 = X(14);
v6 = X(15);
v7 = X(16); %1.0297  Wert aus PSAT (8. Knoten entnommen)
v8 = X(17);
v9 = X(18);


% Vektoren f???r 9 Bus Netzwerk
Pw = [ p1 p2 p3 p4 p5 p6 p7 p8 p9]';
Q =  [ q1 q2 q3 q4 q5 q6 q7 q8 q9]';
V =  [ v1 v2 v3 v4 v5 v6 v7 v8 v9]';
Theta = [theta1 theta2 theta3 theta4 theta5 theta6 theta7 theta8 theta9];

n = length (Pw);

Theta_hk = Theta'*ones(1,n)-ones(n,1)*Theta ; % Thetamatrix mit Winkeldiffernz auf der Diagonale


X0_G = [

V.*( (real(G.Ystart).*cos(Theta_hk) + imag(G.Ystart).*sin(Theta_hk))* V ) - Pw ;...   % Netzwerk

V.*( (real(G.Ystart).*sin(Theta_hk) - imag(G.Ystart).*cos(Theta_hk))* V ) - Q ;... 

];

function X0_gen = A_init_gen (X,P,I)

% Generator 1
delta1 = X(1);
omega1 = X(2);
Pm1    = X(3);
E1     = X(4);

% Generator 2
delta2 = X(5);
omega2 = X(6);
Pm2    = X(7);
E2     = X(8);


% Generator 3
delta3 = X(9);
omega3 = X(10);
Pm3    = X(11);
E3     = X(12);

X0_gen = [

% Generator 1
P.omegaS*(omega1);...
(1/P.M1)*(Pm1-I.Pe1 - (P.D1*omega1));...
I.Pe1 - (E1*I.V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-I.theta1));...
I.Q1  + (E1*I.V1*abs(P.xd1)*sin(angle(P.xd1)+delta1-I.theta1)) - (I.V1^2*abs(P.xd1)*sin(angle(P.xd1)));...

% Generator 2
P.omegaS*(omega2);...
(1/P.M2)*(Pm2-I.Pe2 - (P.D2*omega2));...
I.Pe2 - (E2*I.V2*abs(P.xd2)*cos(angle(P.xd2)+delta2-I.theta2));
I.Q2  + (E2*I.V2*abs(P.xd2)*sin(angle(P.xd2)+delta2-I.theta2)) - (I.V2^2*abs(P.xd2)*sin(angle(P.xd2)));

% Generator 3
% P.omegaS*(omega3-P.omega1);...             %3 DGL =0
P.omegaS*(omega3);...
(1/P.M3)*(Pm3-I.Pe3 - (P.D3*omega3));...
I.Pe3 - (E3*I.V3*abs(P.xd3)*cos(angle(P.xd3)+delta3-I.theta3))
I.Q3  + (E3*I.V3*abs(P.xd3)*sin(angle(P.xd3)+delta3-I.theta3)) - (I.V3^2*abs(P.xd3)*sin(angle(P.xd3)))

];

function out = LPV_3Gen(t,X,xs,I,P)
%% Generatoren und Reglervar.

% Generator 1
delta1=X(1,:);
omega1=X(2,:);


% Generator 2
delta2=X(3,:);
omega2=X(4,:);

% Generator 3
delta3=X(5,:);
omega3=X(6,:);


% Netzwerk
p1 = X(24,:);
p2 = X(25,:);
p3 = X(26,:);
p4 = X(27,:);
p5 = X(28,:);
p6 = X(29,:);
p7 = X(30,:);
p8 = X(31,:);
p9 = X(32,:);

q1 = X(33,:);
q2 = X(34,:);
q3 = X(35,:);
q4 = X(36,:);
q5 = X(37,:);
q6 = X(38,:);
q7 = X(39,:);
q8 = X(40,:);
q9 = X(41,:);

theata2 = X(7,:);
theata3 = X(8,:);
theta4 = X(9,:);
theta5 = X(10,:);
theta6 = X(11,:);
theta7 = X(12,:);
theta8 = X(13,:);
theta9 = X(14,:);

V1  = X(15,:);

V2 = X(16,:);
V3 = X(17,:);
v4 = X(18,:);
v5 = X(19,:);
v6 = X(20,:);
v7 = X(21,:);
v8 = X(22,:);
v9 = X(23,:);
%v9 = X(42,:);


%% Vektoren f??r 9 Bus Netzwerk
V =  [ V1 V2 V3 v4 v5 v6 v7 v8 v9].';

Pw = [ p1 p2 p3 p4 p5 p6 p7 p8 p9].';
Q =  [ q1 q2 q3 q4 q5 q6 q7 q8 q9].';

Theta = [0 theata2 theata3 theta4 theta5 theta6 theta7 theta8 theta9];

n = length (Pw);
Theta_hk = Theta'*ones(1,n)-ones(n,1)*Theta ;

Pe1 = (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-I.theta1));
Pe2 = (I.E2*V2*abs(P.xd2)*cos(angle(P.xd2)+delta2-theata2));
Pe3 = (I.E3*V3*abs(P.xd3)*cos(angle(P.xd3)+delta3-theata3));


out = [
P.omegaS*omega1
(1/P.M1)*(I.Pm1-Pe1 - (P.D1*omega1))

P.omegaS*omega2
(1/P.M2)*(I.Pm2-Pe2 - (P.D2*omega2))

P.omegaS*omega3
(1/P.M3)*(I.Pm3-Pe3 - (P.D3*omega3))

p1 - (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1));
p2 - (I.E2*V2*abs(P.xd2)*cos(angle(P.xd2)+delta2-theata2));
p3 - (I.E3*V3*abs(P.xd3)*cos(angle(P.xd3)+delta3-theata3));
p4 - I.p4 * v4^2 / I.v4^2;
p5 - I.p5 * v5^2 / I.v5^2;
p6 - I.p6 * v6^2 / I.v6^2;
p7 - I.p7 * v7^2 / I.v7^2;
p8 - I.p8 * v8^2 / I.v8^2;
p9 - I.p9 * v9^2 / I.v9^2;

%q1 + (I.E1*V1*abs(P.xd1)*sin(angle(P.xd1)+delta1)) - (V1^2*abs(P.xd1)*sin(angle(P.xd1)));
q2 + (I.E2*V2*abs(P.xd2)*sin(angle(P.xd2)+delta2-theata2)) - (V2^2*abs(P.xd2)*sin(angle(P.xd2)));
q3 + (I.E3*V3*abs(P.xd3)*sin(angle(P.xd3)+delta3-theata3)) - (V3^2*abs(P.xd3)*sin(angle(P.xd3)));
q4 - I.q4 * v4^2 / I.v4^2;
q5 - I.q5 * v5^2 / I.v5^2;
q6 - I.q6 * v6^2 / I.v6^2;
q7 - I.q7 * v7^2 / I.v7^2;
q8 - I.q8 * v8^2 / I.v8^2;
q9 - I.q9 * v9^2 / I.v9^2;

V.*( (real(xs).*cos(Theta_hk) + imag(xs).*sin(Theta_hk))* V ) - Pw; 
V.*( (real(xs).*sin(Theta_hk) - imag(xs).*cos(Theta_hk))* V ) - Q; 
];
