function completed = example_composite_30bus()
% updated: 30-August-2017, AElG
% As realized in the ACC paper
% Compositional transient stability anaylsis of Power Systems
clear
clc

global P I
load busIEEE30_Y
[P,~,I]   = GetStarted(Ybus);
dim_x = 12;
options = SetOptions(2);


y_gen = [ I.delta1 I.omega1 I.delta2 I.omega2 I.delta5 I.omega5 I.delta8 I.omega8 I.delta11 I.omega11 I.delta13 I.omega13];

        
y_grid = [  I.theta2 I.theta3 I.theta4 I.theta5 I.theta6 I.theta7 I.theta8 I.theta9 I.theta10 I.theta11 I.theta12 I.theta13 I.theta14 ...
            I.theta15 I.theta16 I.theta17 I.theta18 I.theta19 I.theta20 I.theta21 I.theta22 I.theta23 I.theta24 I.theta25 I.theta26 I.theta27 I.theta28 I.theta29 I.theta30 ...
            I.V1 I.V2 I.v3 I.v4 I.V5 I.v6 I.v7 I.V8 I.v9 I.v10 I.V11 I.v12 I.V13 I.v14 I.v15 ...
            I.v16 I.v17 I.v18 I.v19 I.v20 I.v21 I.v22 I.v23 I.v24 I.v25 I.v26 I.v27 I.v28 I.v29 I.v30 ...
            I.Pe1 I.Pe2 I.p3 I.p4 I.Pe5 I.p6 I.p7 I.Pe8 I.p9 I.p10 I.Pe11 I.p12 I.Pe13 I.p14 I.p15 ...
            I.p16 I.p17 I.p18 I.p19 I.p20 I.p21 I.p22 I.p23 I.p24 I.p25 I.p26 I.p27 I.p28 I.p29 I.p30 ...
            I.Q1 I.Q2 I.q3 I.q4 I.Q5 I.q6 I.q7 I.Q8 I.q9 I.q10 I.Q11 I.q12 I.Q13 I.q14 I.q15 ...
            I.q16 I.q17 I.q18 I.q19 I.q20 I.q21 I.q22 I.q23 I.q24 I.q25 I.q26 I.q27 I.q28 I.q29 I.q30];
            
        
y0 = [y_gen...                  % Zust??nde der 3 Generatoren              % Algebr. Variablen der Gen.
      y_grid...                 % Netzvariablen
      ]';
  
%solveDAE(y_gen,y0,P,I);
  
%% Assign parameters for reachability analysis ----------------------------
x0{1} = y_gen(1,1:2);
x0{2} = y_gen(1,3:4);
x0{3} = y_gen(1,5:6);
x0{4} = y_gen(1,7:8);
x0{5} = y_gen(1,9:10);
x0{6} = y_gen(1,11:12);

g_R0  = 1e-3*eye(dim_x,dim_x);

i = 1;
%--------------------------------------------------------------------------
%% Load classes
powerDyn_Gen1  = nonlinearSys(2,1,@Gen1_30,options);
powerDyn_Gen2  = nonlinearSys(2,2,@Gen2_30,options);
powerDyn_Gen5  = nonlinearSys(2,2,@Gen5_30,options);
powerDyn_Gen8  = nonlinearSys(2,2,@Gen8_30,options);
powerDyn_Gen11 = nonlinearSys(2,2,@Gen11_30,options);
powerDyn_Gen13 = nonlinearSys(2,2,@Gen13_30,options);


t                = 0;

Rnext1 = zonotope([x0{1}.',g_R0(1:2,:)]);
Z01    = Rnext1;
Rnext2 = zonotope([x0{2}.',g_R0(3:4,:)]);
Z02    = Rnext2;
Rnext5 = zonotope([x0{3}.',g_R0(5:6,:)]); % x0
Z03    = Rnext5;
Rnext8 = zonotope([x0{4}.',g_R0(7:8,:)]); % x0
Z04    = Rnext8;
Rnext11 = zonotope([x0{5}.',g_R0(9:10,:)]); % x0
Z05    = Rnext11;
Rnext13 = zonotope([x0{6}.',g_R0(11:12,:)]); % x0
Z06    = Rnext13;

options1 = options;
options2 = options;
options5 = options;
options8 = options;
options11 = options;
options13 = options;
options_fsolve = optimoptions('fsolve','Display','off');
tic
if options.reachflag
%% normal operation
    while t < options.tFinal
        
        % Check fault conditions
        if t > 0.01 && t < 0.02
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

        y_gen5   = center(Rnext5);
        G.delta5 = y_gen5(1);
        
        y_gen8   = center(Rnext8);
        G.delta8 = y_gen8(1);
        
        y_gen11   = center(Rnext11);
        G.delta11 = y_gen11(1);
        
        y_gen13   = center(Rnext13);
        G.delta13 = y_gen13(1);

        init_grid = @(X)PowerGrid(X,G,P,I);
        X0        = fsolve(init_grid, [zeros(29,1);ones(90,1)],options_fsolve );
        
        
        % Gen1
        options1.R0        = Rnext1; % x0
        options1.uTrans    = X0(30,1);
        options1.U         = zonotope([0 0.0085]);

        [Rcont1,Rnext1,options1] = contReach(powerDyn_Gen1,options1);

        % Gen2
        options2.R0        = Rnext2; % x0
        options2.uTrans    = [X0(1,1);X0(31,1)];
        options2.U         = zonotope([zeros(2,1) [0.00025 0;0 0.00015]]);
        
        [Rcont2,Rnext2,options2] = contReach(powerDyn_Gen2,options2);

        % Gen5
        options5.R0        = Rnext5; % x0
        options5.uTrans    = [X0(4,1);X0(34,1)];
        options5.U         = zonotope([zeros(2,1) [0.00035 0;0 0.0025]]);

        [Rcont5,Rnext5,options5] = contReach(powerDyn_Gen5,options5);
        
        % Gen8
        options8.R0        = Rnext8; % x0
        options8.uTrans    = [X0(7,1);X0(37,1)];
        options8.U         = zonotope([zeros(2,1) [0.00035 0;0 0.0025]]);
       
        [Rcont8,Rnext8,options8] = contReach(powerDyn_Gen8,options8);
        
        % Gen11
        options11.R0        = Rnext11; % x0
        options11.uTrans    = [X0(10,1);X0(40,1)];
        options11.U         = zonotope([zeros(2,1) [0.0005 0;0 0.004]]);
        

        [Rcont11,Rnext11,options11] = contReach(powerDyn_Gen11,options11);
        
        % Gen13
        options13.R0        = Rnext13; % x0
        options13.uTrans    = [X0(12,1);X0(42,1)];
        options13.U         = zonotope([zeros(2,1) [0.00055 0;0 0.0035]]);
        
        [Rcont13,Rnext13,options13] = contReach(powerDyn_Gen13,options13);
        
        Rcont{i}    = cartesianProduct(Rcont1,Rcont2);
        Rcont{i}    = cartesianProduct(Rcont{i},Rcont5);
        Rcont{i}    = cartesianProduct(Rcont{i},Rcont8);
        Rcont{i}    = cartesianProduct(Rcont{i},Rcont11);
        Rcont{i}    = cartesianProduct(Rcont{i},Rcont13);
        
        subplot(2,3,1)
        hold on
        plotFilled(Rcont{i},[1 2],[.8 .8 .8],'EdgeColor','none');
        subplot(2,3,2)
        hold on
        plotFilled(Rcont{i},[3 4],[.8 .8 .8],'EdgeColor','none');
        subplot(2,3,3)
        hold on
        plotFilled(Rcont{i},[5 6],[.8 .8 .8],'EdgeColor','none');
        subplot(2,3,4)
        hold on
        plotFilled(Rcont{i},[7 8],[.8 .8 .8],'EdgeColor','none');
        subplot(2,3,5)
        hold on
        plotFilled(Rcont{i},[9 10],[.8 .8 .8],'EdgeColor','none');
        subplot(2,3,6)
        hold on
        plotFilled(Rcont{i},[11 12],[.8 .8 .8],'EdgeColor','none');
        
         
        t = t + options.timeStep;
        t
         
       
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
        y_gen4=randPointExtreme(Z04); % At the borders
        y_gen5=randPointExtreme(Z05); % At the borders
        y_gen6=randPointExtreme(Z06); % At the borders
    else
        y_gen1=randPoint(Z01); % Inside the Zonotope
        y_gen2=randPoint(Z02); % Inside the Zonotope
        y_gen3=randPoint(Z03); % Inside the Zonotope
        y_gen4=randPoint(Z04); % Inside the Zonotope
        y_gen5=randPoint(Z05); % Inside the Zonotope
        y_gen6=randPoint(Z06); % Inside the Zonotope
           
    end

    
y_grid = [  I.theta2 I.theta3 I.theta4 I.theta5 I.theta6 I.theta7 I.theta8 I.theta9 I.theta10 I.theta11 I.theta12 I.theta13 I.theta14 ...
            I.theta15 I.theta16 I.theta17 I.theta18 I.theta19 I.theta20 I.theta21 I.theta22 I.theta23 I.theta24 I.theta25 I.theta26 I.theta27 I.theta28 I.theta29 I.theta30 ...
            I.V1 I.V2 I.v3 I.v4 I.V5 I.v6 I.v7 I.V8 I.v9 I.v10 I.V11 I.v12 I.V13 I.v14 I.v15 ...
            I.v16 I.v17 I.v18 I.v19 I.v20 I.v21 I.v22 I.v23 I.v24 I.v25 I.v26 I.v27 I.v28 I.v29 I.v30 ...
            I.Pe1 I.Pe2 I.p3 I.p4 I.Pe5 I.p6 I.p7 I.Pe8 I.p9 I.p10 I.Pe11 I.p12 I.Pe13 I.p14 I.p15 ...
            I.p16 I.p17 I.p18 I.p19 I.p20 I.p21 I.p22 I.p23 I.p24 I.p25 I.p26 I.p27 I.p28 I.p29 I.p30 ...
            I.Q1 I.Q2 I.q3 I.q4 I.Q5 I.q6 I.q7 I.Q8 I.q9 I.q10 I.Q11 I.q12 I.Q13 I.q14 I.q15 ...
            I.q16 I.q17 I.q18 I.q19 I.q20 I.q21 I.q22 I.q23 I.q24 I.q25 I.q26 I.q27 I.q28 I.q29 I.q30];
            
y_gen = [y_gen1;y_gen2;y_gen3;y_gen4;y_gen5;y_gen6].'  ;  
y0 = [ y_gen y_grid].';

solveDAE(y_gen,y0,P,I)
end

subplot(2,3,1)
hold on
plotFilled(Z01,[1 2],[1 1 1],'EdgeColor',[.2 .2 .2]); 
subplot(2,3,2)
hold on
plotFilled(Z02,[1 2],[1 1 1],'EdgeColor',[.2 .2 .2]); 
subplot(2,3,3)
hold on
plotFilled(Z03,[1 2],[1 1 1],'EdgeColor',[.2 .2 .2]); 
subplot(2,3,4)
hold on
plotFilled(Z04,[1 2],[1 1 1],'EdgeColor',[.2 .2 .2]); 
subplot(2,3,5)
hold on
plotFilled(Z05,[1 2],[1 1 1],'EdgeColor',[.2 .2 .2]); 
subplot(2,3,6)
hold on
plotFilled(Z06,[1 2],[1 1 1],'EdgeColor',[.2 .2 .2]); 

%example completed
completed = 1;

function X0_G = PowerGrid(x,G,P,I)

xs       = G.xs;
delta1   = G.delta1;
delta2   = G.delta2;
delta5   = G.delta5;
delta8   = G.delta8;
delta11  = G.delta11;
delta13  = G.delta13;

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
theata10 = x(9);
theata11 = x(10);
theata12 = x(11);
theata13 = x(12);
theata14 = x(13);
theata15 = x(14);
theata16 = x(15);
theata17 = x(16);
theata18 = x(17);
theata19 = x(18);
theata20 = x(19);
theata21 = x(20);
theata22 = x(21);
theata23 = x(22);
theata24 = x(23);
theata25 = x(24);
theata26 = x(25);
theata27 = x(26);
theata28 = x(27);
theata29 = x(28);
theata30 = x(29);

% Bus voltage levels
V1  = x(30);
V2  = x(31);
v3  = x(32);
v4  = x(33);
V5  = x(34);
v6  = x(35);
v7  = x(36);
V8  = x(37);
v9  = x(38);
v10 = x(39);
V11 = x(40);
v12 = x(41);
V13 = x(42);
v14 = x(43);
v15 = x(44);
v16 = x(45);
v17 = x(46);
v18 = x(47);
v19 = x(48);
v20 = x(49);
v21 = x(50);
v22 = x(51);
v23 = x(52);
v24 = x(53);
v25 = x(54);
v26 = x(55);
v27 = x(56);
v28 = x(57);
v29 = x(58);
v30 = x(59);


p1  = x(60);
p2  = x(61);
p3  = x(62);
p4  = x(63);
p5  = x(64); 
p6  = x(65);
p7  = x(66);
p8  = x(67);
p9  = x(68);
p10 = x(69);
p11 = x(70);
p12 = x(71);
p13 = x(72);
p14 = x(73);
p15 = x(74);
p16 = x(75);
p17 = x(76);
p18 = x(77);
p19 = x(78);
p20 = x(79);
p21 = x(80);
p22 = x(81);
p23 = x(82);
p24 = x(83);
p25 = x(84);
p26 = x(85);
p27 = x(86);
p28 = x(87);
p29 = x(88);
p30 = x(89);

q1  = x(90);
q2  = x(91);
q3  = x(92);
q4  = x(93);
q5  = x(94);
q6  = x(95);
q7  = x(96);
q8  = x(97);
q9  = x(98);
q10 = x(99);
q11 = x(100);
q12 = x(101);
q13 = x(102);
q14 = x(103);
q15 = x(104);
q16 = x(105);
q17 = x(106);
q18 = x(107);
q19 = x(108);
q20 = x(109);
q21 = x(110);
q22 = x(111);
q23 = x(112);
q24 = x(113);
q25 = x(114);
q26 = x(115);
q27 = x(116);
q28 = x(117);
q29 = x(118);
q30 = x(119);

V =  [ V1 V2 v3 v4 V5 v6 v7 V8 v9 v10 V11 v12 V13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 ].';

Pw = [ p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30].';
Q =  [ q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16 q17 q18 q19 q20 q21 q22 q23 q24 q25 q26 q27 q28 q29 q30].';

Theta = [0 theata2 theata3 theata4 theata5 theata6 theata7 theata8 theata9 theata10 theata11 theata12 theata13 theata14 theata15 ...
           theata16 theata17 theata18 theata19 theata20 theata21 theata22 theata23 theata24 theata25 theata26 theata27 theata28 theata29 theata30];

n        = length (V);
Theta_hk = Theta.'*ones(1,n)-ones(n,1)*Theta ;

X0_G = [ p1 - (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-theata1))
         p2 - (I.E2*V2*abs(P.xd2)*cos(angle(P.xd2)+delta2-theata2));
         p3 - I.p3 * v3^2 / I.v3^2;
p4 - I.p4 * v4^2 / I.v4^2;
p5 - (I.E5*V5*abs(P.xd5)*cos(angle(P.xd5)+delta5-theata5));
p6 - I.p6 * v6^2 / I.v6^2;
p7 - I.p7 * v7^2 / I.v7^2;
p8 - (I.E8*V8*abs(P.xd8)*cos(angle(P.xd8)+delta8-theata8));
p9 - I.p9 * v9^2 / I.v9^2;
p10 - I.p10 * v10^2 / I.v10^2;
p11 - (I.E11*V11*abs(P.xd11)*cos(angle(P.xd11)+delta11-theata11));
p12 - I.p12 * v12^2 / I.v12^2;
p13 - (I.E13*V13*abs(P.xd13)*cos(angle(P.xd13)+delta13-theata13));
p14 - I.p14 * v14^2 / I.v14^2;
p15 - I.p15 * v15^2 / I.v15^2;
p16 - I.p16 * v16^2 / I.v16^2;
p17 - I.p17 * v17^2 / I.v17^2;
p18 - I.p18 * v18^2 / I.v18^2;
p19 - I.p19 * v19^2 / I.v19^2;
p20 - I.p20 * v20^2 / I.v20^2;
p21 - I.p21 * v21^2 / I.v21^2;
p22 - I.p22 * v22^2 / I.v22^2;
p23 - I.p23 * v23^2 / I.v23^2;
p24 - I.p24 * v24^2 / I.v24^2;
p25 - I.p25 * v25^2 / I.v25^2;
p26 - I.p26 * v26^2 / I.v26^2;
p27 - I.p27 * v27^2 / I.v27^2;
p28 - I.p28 * v28^2 / I.v28^2;
p29 - I.p29 * v29^2 / I.v29^2;
p30 - I.p30 * v30^2 / I.v30^2;

%dx(10)  = q1 - (I.E1*V1*abs(P.xd1)*sin(angle(P.xd1)+delta1-theata1)) + (V1^2*abs(P.xd1)*sin(angle(P.xd1)));
q2 + (I.E2*V2*abs(P.xd2)*sin(angle(P.xd2)+delta2-theata2)) - (V2^2*abs(P.xd2)*sin(angle(P.xd2)));
q3 - I.q3 * v3^2 / I.v3^2;
q4 - I.q4 * v4^2 / I.v4^2;
q5 + (I.E5*V5*abs(P.xd5)*sin(angle(P.xd5)+delta5-theata5)) - (V5^2*abs(P.xd5)*sin(angle(P.xd5)));
q6 - I.q6 * v6^2 / I.v6^2;
q7 - I.q7 * v7^2 / I.v7^2;
q8 + (I.E8*V8*abs(P.xd8)*sin(angle(P.xd8)+delta8-theata8)) - (V8^2*abs(P.xd8)*sin(angle(P.xd8)));
q9 - I.q9 * v9^2 / I.v9^2;
q10 - I.q10 * v10^2 / I.v10^2;
q11 + (I.E11*V11*abs(P.xd11)*sin(angle(P.xd11)+delta11-theata11)) - (V11^2*abs(P.xd11)*sin(angle(P.xd11)));
q12 - I.q12 * v12^2 / I.v12^2;
q13 + (I.E13*V13*abs(P.xd13)*sin(angle(P.xd13)+delta13-theata13)) - (V13^2*abs(P.xd13)*sin(angle(P.xd13)));
q14 - I.q14 * v14^2 / I.v14^2;
q15 - I.q15 * v15^2 / I.v15^2;
q16 - I.q16 * v16^2 / I.v16^2;
q17 - I.q17 * v17^2 / I.v17^2;
q18 - I.q18 * v18^2 / I.v18^2;
q19 - I.q19 * v19^2 / I.v19^2;
q20 - I.q20 * v20^2 / I.v20^2;
q21 - I.q21 * v21^2 / I.v21^2;
q22 - I.q22 * v22^2 / I.v22^2;
q23 - I.q23 * v23^2 / I.v23^2;
q24 - I.q24 * v24^2 / I.v24^2;
q25 - I.q25 * v25^2 / I.v25^2;
q26 - I.q26 * v26^2 / I.v26^2;
q27 - I.q27 * v27^2 / I.v27^2;
q28 - I.q28 * v28^2 / I.v28^2;
q29 - I.q29 * v29^2 / I.v29^2;
q30 - I.q30 * v30^2 / I.v30^2;

V.*( (real(xs).*cos(Theta_hk) + imag(xs).*sin(Theta_hk))* V ) - Pw; 
V.*( (real(xs).*sin(Theta_hk) - imag(xs).*cos(Theta_hk))* V ) - Q; 
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
t3 = 0.02;
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
f1r = @(t,x)IEEE30(t,x,xs,I,P);
xs = P.xs1;
f2r = @(t,x)IEEE30(t,x,xs,I,P);
xs = P.xs;
f3r = @(t,x)IEEE30(t,x,xs,I,P);

[~,y1r] = ode15s(f1r,tspan1,y0,options_fsolve);
[~,y2r] = ode15s(f2r,tspan2,y1r(end,:),options_fsolve);
[~,y3r] = ode15s(f3r,tspan3,y2r(end,:),options_fsolve);

%Tr=[t1r;t2r;t3r];

Yr=[y1r;y2r;y3r];

%Ef=Yr(:,[61, 62,63]);


figure(1)
subplot(2,3,1)
hold on
plot(Yr(:,1), Yr(:,2),'k')
subplot(2,3,2)
hold on
plot(Yr(:,3), Yr(:,4),'k')
subplot(2,3,3)
hold on
plot(Yr(:,5),Yr(:,6),'k')
subplot(2,3,4)
hold on
plot(Yr(:,7), Yr(:,8),'k')
subplot(2,3,5)
hold on
plot(Yr(:,9), Yr(:,10),'k')
subplot(2,3,6)
hold on
plot(Yr(:,11), Yr(:,12),'k')

function [options] = SetOptions(dim_x)
%global P

options.tStart   = 0;   % Start time
options.tFinal   = 0.5;   % Final time 
options.timeStep = 0.001;

options.Type = 'Power System';
options.reachflag = 1;

% Set options -------------------------------------------------------------


%% Initial set for the input variables ------------------------------------
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

function [P,G,I] = GetStarted(Ybus30)
%% Assign fixed parameters of the synchronous generator
% 1. Generator
P.M1         = 1/(15*pi);   	% Tr???gheitskoeffizient
P.D1         = 1;           % D???mpfungskoeffizient
P.xd1        = 1/(1i*0.2);  	% d-Achsen transiente Reaktanz     


% 2. Generator
P.M2         = 1/(15*pi);   	 % Tr???gheitskoeffizient
P.D2         = 1;        % D???mpfungskoeffizient
P.xd2        = 1/(1i*0.2);  	% d-Achsen transiente Reaktanz

% 3. Generator
P.M5         = 1/(15*pi);   	% Tr???gheitskoeffizient
P.D5         = 1;        % D???mpfungskoeffizient
P.xd5        = 1/(1i*0.2);  	% d-Achsen transiente Reaktanz

% 4. Generator
P.M8         = 1/(15*pi);   	% Tr???gheitskoeffizient
P.D8         = 1;        % D???mpfungskoeffizient
P.xd8        = 1/(1i*0.2);  	% d-Achsen transiente Reaktanz

% 5. Generator
P.M11         = 1/(15*pi);   	% Tr???gheitskoeffizient
P.D11         = 1;        % D???mpfungskoeffizient
P.xd11        = 1/(1i*0.2);  	% d-Achsen transiente Reaktanz

% 6. Generator
P.M13         = 1/(15*pi);   	% Tr???gheitskoeffizient
P.D13         = 1;        % D???mpfungskoeffizient
P.xd13        = 1/(1i*0.2);  	% d-Achsen transiente Reaktanz

P.omegaS     = 2*pi*50;

%% Network parameters

% Definition of the ADMITTANCE MATRIX --- Aufstellen der Admittanzmatrix f???r den Power Flow
% Original-System as describes by Anderson and FOuad --- Originalsystem
G.Y         = full(Ybus30);

%% Admittance-matrizes for different faults
% Admittance change of line from Bus 5 to Bus 7: half of the admittances
G.Yf              = G.Y;
G.Yf(2,4)         = G.Yf(2,6)/1.1;
G.Yf(4,2)         = G.Yf(6,2)/1.1;


xs = G.Yf;        % Originale Admittanzmatrix 
G.Ystart=xs;
xs1 = G.Y;        % Admittanzmatrix mit Fehler G.Yf=nur Halbierung der Leitung, G.Yf1=Schlimm

P.xs  = xs;  % Normal operation
P.xs1 = xs1; % Fault


%% Power Flow: Startvorgaben
G.theta1    = 0;             % Winkel am Slack Bus

G.v1         = 1.0600;         % Spannung am Slack Bus (Referenzknoten 1)
G.v2         = 1.043;         % Spannung an Gen. 2
G.v5         = 1.010;         % Spannung an Gen. 3
G.v8         = 1.010;
G.v11        = 1.082;
G.v13        = 1.071;

G.p2         = 0.5;          % Wirkleistung an Gen. 2
G.p5         = 0.5;          % Wirkleistung an Gen. 2
G.p8         = 0.5;          % Wirkleistung an Gen. 2
G.p11        = 0.5;          % Wirkleistung an Gen. 2
G.p13        = 0.5;          % Wirkleistung an Gen. 2
    
    
G.p3  = -0.0240;
G.p4  = -0.0760;
G.p6  =  0;
G.p7  = -0.2280;
G.p9  =  0;
G.p10 = -0.0580;
G.p12 = -0.1120;
G.p14 = -0.0620;
G.p15 = -0.0820;
G.p16 = -0.0350;
G.p17 = -0.0900;
G.p18 = -0.0320;
G.p19 = -0.0950;
G.p20 = -0.0220;
G.p21 = -0.1750;
G.p22 = 0;
G.p23 = -0.0320;
G.p24 = -0.0870;
G.p25 = 0;
G.p26 = -0.0350;
G.p27 = 0;
G.p28 = 0;
G.p29 = -0.0240;
G.p30 = -0.1060;



G.q3   = -0.0120;
G.q4   = -0.0160;
G.q6   = 0;
G.q7   = -0.1090;
G.q9   = 0;
G.q10  = -0.0200;
G.q12  = -0.0750;
G.q14  = -0.0160;
G.q15  = -0.0250;
G.q16  = -0.0180;
G.q17  = -0.0580;
G.q18  = -0.0090;
G.q19  = -0.0340;
G.q20  = -0.0070;
G.q21  = -0.1120;
G.q22  = -0;
G.q23  = -0.0160;
G.q24  = -0.0670;
G.q25  = -0;
G.q26  = -0.0230;
G.q27  = 0;
G.q28  = 0;
G.q29  = -0.0090;
G.q30  = -0.0190;


%Solution of the Power Flows Problem
options = optimoptions('fsolve','Display','off');

init_grid = @(X)A_init_grid_30k(X,G);
X0 = fsolve(init_grid, ones(60,1),options);


% Initial Values of the grid after Power Flow
I.Pe1  = X0(1);
I.Pe2  = G.p2;
I.p3   = G.p3;
I.p4   = G.p4;
I.Pe5  = G.p5;
I.p6   = G.p6;
I.p7   = G.p7; 
I.Pe8  = G.p8;
I.p9   = G.p9;
I.p10  = G.p10;
I.Pe11 = G.p11;
I.p12  = G.p12;
I.Pe13 = G.p13;
I.p14  = G.p14;
I.p15  = G.p15;
I.p16  = G.p16;
I.p17  = G.p17;
I.p18  = G.p18;
I.p19  = G.p19;
I.p20  = G.p20;
I.p21  = G.p21;
I.p22  = G.p22;
I.p23  = G.p23;
I.p24  = G.p24;
I.p25  = G.p25;
I.p26  = G.p26;
I.p27  = G.p27;
I.p28  = G.p28;
I.p29  = G.p29;
I.p30  = G.p30;



I.Q1   = X0(2);
I.Q2   = X0(3);
I.q3   = G.q3;
I.q4   = G.q4;
I.Q5   = X0(4);
I.q6   = G.q6;
I.q7   = G.q7;
I.Q8   = X0(5);
I.q9   = G.q9;
I.q10  = G.q10;
I.Q11  = X0(6);
I.q12  = G.q12;
I.Q13  = X0(7);
I.q14  = G.q14;
I.q15  = G.q15;
I.q16  = G.q16;
I.q17  = G.q17;
I.q18  = G.q18;
I.q19  = G.q19;
I.q20  = G.q20;
I.q21  = G.q21;
I.q22  = G.q22;
I.q23  = G.q23;
I.q24  = G.q24;
I.q25  = G.q25;
I.q26  = G.q26;
I.q27  = G.q27;
I.q28  = G.q28;
I.q29  = G.q29;
I.q30  = G.q30;


I.theta1  = G.theta1;
I.theta2  = X0(8);
I.theta3  = X0(9);
I.theta4  = X0(10);
I.theta5  = X0(11);
I.theta6  = X0(12);
I.theta7  = X0(13); 
I.theta8  = X0(14);
I.theta9  = X0(15);
I.theta10 = X0(16);
I.theta11 = X0(17);
I.theta12 = X0(18);
I.theta13 = X0(19);
I.theta14 = X0(20);
I.theta15 = X0(21);
I.theta16 = X0(22);
I.theta17 = X0(23);
I.theta18 = X0(24);
I.theta19 = X0(25);
I.theta20 = X0(26);
I.theta21 = X0(27);
I.theta22 = X0(28);
I.theta23 = X0(29);
I.theta24 = X0(30);
I.theta25 = X0(31);
I.theta26 = X0(32);
I.theta27 = X0(33);
I.theta28 = X0(34);
I.theta29 = X0(35);
I.theta30 = X0(36);


I.V1  = G.v1;
I.V2  = G.v2;
I.v3  = X0(37);
I.v4  = X0(38);
I.V5  = G.v5;
I.v6  = X0(39);
I.v7  = X0(40); 
I.V8  = G.v8;
I.v9  = X0(41);
I.v10 = X0(42);
I.V11 = G.v11;
I.v12 = X0(43);
I.V13 = G.v13;
I.v14 = X0(44);
I.v15 = X0(45);
I.v16 = X0(46);
I.v17 = X0(47);
I.v18 = X0(48);
I.v19 = X0(49);
I.v20 = X0(50);
I.v21 = X0(51);
I.v22 = X0(52);
I.v23 = X0(53);
I.v24 = X0(54);
I.v25 = X0(55);
I.v26 = X0(56);
I.v27 = X0(57);
I.v28 = X0(58);
I.v29 = X0(59);
I.v30 = X0(60);


alg = @(X)A_init_gen(X,P,I);
X0_gen = fsolve(alg, ones(24,1),options);

% Generator 1
I.delta1 = X0_gen(1);
I.omega1 = X0_gen(2);
I.Pm1    = X0_gen(3);
I.E1     = X0_gen(4);

% Generator 2
I.delta2 = X0_gen(5);
I.omega2 = X0_gen(6);
I.Pm2    = X0_gen(7);
I.E2     = X0_gen(8);

% Generator 5
I.delta5 = X0_gen(9);
I.omega5 = X0_gen(10);
I.Pm5    = X0_gen(11);
I.E5     = X0_gen(12);

% Generator 8
I.delta8 = X0_gen(13);
I.omega8 = X0_gen(14);
I.Pm8    = X0_gen(15);
I.E8     = X0_gen(16);

% Generator 1
I.delta11 = X0_gen(17);
I.omega11 = X0_gen(18);
I.Pm11    = X0_gen(19);
I.E11     = X0_gen(20);

% Generator 8
I.delta13 = X0_gen(21);
I.omega13 = X0_gen(22);
I.Pm13    = X0_gen(23);
I.E13     = X0_gen(24);

function X0_G = A_init_grid_30k (X,G)

% Parameter f???r das Netz
p1  = X(1);
p2  = G.p2;
p3  = G.p3;
p4  = G.p4;
p5  = G.p5;
p6  = G.p6;
p7  = G.p7;  % X(13)
p8  = G.p8;
p9  = G.p9;
p10 = G.p10;
p11 = G.p11;
p12 = G.p12;
p13 = G.p13;
p14 = G.p14;
p15 = G.p15;
p16 = G.p16;
p17 = G.p17;
p18 = G.p18;
p19 = G.p19;
p20 = G.p20;
p21 = G.p21;
p22 = G.p22;
p23 = G.p23;
p24 = G.p24;
p25 = G.p25;
p26 = G.p26;
p27 = G.p27;
p28 = G.p28;
p29 = G.p29;
p30 = G.p30;

q1  = X(2);
q2  = X(3);
q3  = G.q3;
q4  = G.q4;
q5  = X(4);
q6  = G.q6;
q7  = G.q7;  % ;X(14)
q8  = X(5);
q9  = G.q9;
q10 = G.q10;
q11 = X(6);
q12 = G.q12;
q13 = X(7);
q14 = G.q14;
q15 = G.q15;
q16 = G.q16;
q17 = G.q17;
q18 = G.q18;
q19 = G.q19;
q20 = G.q20;
q21 = G.q21;
q22 = G.q22;
q23 = G.q23;
q24 = G.q24;
q25 = G.q25;
q26 = G.q26;
q27 = G.q27;
q28 = G.q28;
q29 = G.q29;
q30 = G.q30;

theta1  = G.theta1;
theta2  = X(8);
theta3  = X(9);
theta4  = X(10);
theta5  = X(11);
theta6  = X(12);
theta7  = X(13); %0.21431 Wert aus PSAT (8. Knoten entnommen)
theta8  = X(14);
theta9  = X(15);
theta10 = X(16);
theta11 = X(17);
theta12 = X(18);
theta13 = X(19);
theta14 = X(20);
theta15 = X(21);
theta16 = X(22);
theta17 = X(23);
theta18 = X(24);
theta19 = X(25);
theta20 = X(26);
theta21 = X(27);
theta22 = X(28);
theta23 = X(29);
theta24 = X(30);
theta25 = X(31);
theta26 = X(32);
theta27 = X(33);
theta28 = X(34);
theta29 = X(35);
theta30 = X(36);

v1 = G.v1;
v2 = G.v2;
v3 = X(37);
v4 = X(38);
v5 = G.v5;
v6 = X(39);
v7 = X(40); %1.0297  Wert aus PSAT (8. Knoten entnommen)
v8 = G.v8;
v9 = X(41);
v10 = X(42);
v11 = G.v11;
v12 = X(43);
v13 = G.v13; 
v14 = X(44);
v15 = X(45);
v16 = X(46);
v17 = X(47);
v18 = X(48);
v19 = X(49);
v20 = X(50);
v21 = X(51);
v22 = X(52);
v23 = X(53);
v24 = X(54);
v25 = X(55);
v26 = X(56);
v27 = X(57);
v28 = X(58);
v29 = X(59);
v30 = X(60);




% Vektoren f???r 9 Bus Netzwerk
Pw = [ p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30]';
Q =  [ q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16 q17 q18 q19 q20 q21 q22 q23 q24 q25 q26 q27 q28 q29 q30]';
V =  [ v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30]';
Theta = [theta1 theta2 theta3 theta4 theta5 theta6 theta7 theta8 theta9 theta10 theta11 theta12 theta13 theta14 theta15...
        theta16 theta17 theta18 theta19 theta20 theta21 theta22 theta23 theta24 theta25 theta26 theta27 theta28 theta29 theta30];

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
delta5 = X(9);
omega5 = X(10);
Pm5    = X(11);
E5     = X(12);


% Generator 4
delta8 = X(13);
omega8 = X(14);
Pm8    = X(15);
E8     = X(16);

% Generator 5
delta11 = X(17);
omega11 = X(18);
Pm11    = X(19);
E11     = X(20);

% Generator 6
delta13 = X(21);
omega13 = X(22);
Pm13    = X(23);
E13     = X(24);

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
P.omegaS*(omega5);...
(1/P.M5)*(Pm5-I.Pe5 - (P.D5*omega5));...
I.Pe5 - (E5*I.V5*abs(P.xd5)*cos(angle(P.xd5)+delta5-I.theta5))
I.Q5  + (E5*I.V5*abs(P.xd5)*sin(angle(P.xd5)+delta5-I.theta5)) - (I.V5^2*abs(P.xd5)*sin(angle(P.xd5)))

% Generator 4
% P.omegaS*(omega3-P.omega1);...             %3 DGL =0
P.omegaS*(omega8);...
(1/P.M8)*(Pm8-I.Pe8 - (P.D8*omega8));...
I.Pe8 - (E8*I.V8*abs(P.xd8)*cos(angle(P.xd8)+delta8-I.theta8))
I.Q8  + (E8*I.V8*abs(P.xd8)*sin(angle(P.xd8)+delta8-I.theta8)) - (I.V8^2*abs(P.xd8)*sin(angle(P.xd8)))

% Generator 5
% P.omegaS*(omega3-P.omega1);...             %3 DGL =0
P.omegaS*(omega11);...
(1/P.M11)*(Pm11-I.Pe11 - (P.D11*omega11));...
I.Pe11 - (E11*I.V11*abs(P.xd11)*cos(angle(P.xd11)+delta11-I.theta11))
I.Q11  + (E11*I.V11*abs(P.xd11)*sin(angle(P.xd11)+delta11-I.theta11)) - (I.V11^2*abs(P.xd11)*sin(angle(P.xd11)))

% Generator 6
% P.omegaS*(omega3-P.omega1);...             %3 DGL =0
P.omegaS*(omega13);...
(1/P.M13)*(Pm13-I.Pe13 - (P.D13*omega13));...
I.Pe13 - (E13*I.V13*abs(P.xd13)*cos(angle(P.xd13)+delta13-I.theta13))
I.Q13  + (E13*I.V13*abs(P.xd13)*sin(angle(P.xd13)+delta13-I.theta13)) - (I.V13^2*abs(P.xd13)*sin(angle(P.xd13)))

];

function out = IEEE30(t,X,xs,I,P)
%% Generatoren und Reglervar.
%global P I
%xs = P.xs;

%if iChange == 2
    %xs = P.xs1;
%end

% Generator 1
delta1=X(1,:);
omega1=X(2,:);


% Generator 2
delta2=X(3,:);
omega2=X(4,:);

% Generator 3
delta5=X(5,:);
omega5=X(6,:);

% Generator 6
delta8=X(7,:);
omega8=X(8,:);

% Generator 8
delta11=X(9,:);
omega11=X(10,:);

% Generator 8
delta13=X(11,:);
omega13=X(12,:);

theata2  = X(13,:);
theta3   = X(14,:);
theta4   = X(15,:);
theata5  = X(16,:);
theta6   = X(17,:);
theta7   = X(18,:);
theata8  = X(19,:);
theta9   = X(20,:);
theta10  = X(21,:);
theata11 = X(22,:);
theta12  = X(23,:);
theata13 = X(24,:);
theta14  = X(25,:);
theta15  = X(26,:);
theta16  = X(27,:);
theta17  = X(28,:);
theta18  = X(29,:);
theta19  = X(30,:);
theta20  = X(31,:);
theta21  = X(32,:);
theta22  = X(33,:);
theta23  = X(34,:);
theta24  = X(35,:);
theta25  = X(36,:);
theta26  = X(37,:);
theta27  = X(38,:);
theta28  = X(39,:);
theta29  = X(40,:);
theta30  = X(41,:);


V1  = X(42,:);

V2  = X(43,:);
v3  = X(44,:);
v4  = X(45,:);
V5  = X(46,:);
v6  = X(47,:);
v7  = X(48,:);
V8  = X(49,:);
v9  = X(50,:);
v10 = X(51,:);
V11 = X(52,:);
v12 = X(53,:);
V13 = X(54,:);
v14 = X(55,:);
v15 = X(56,:);
v16 = X(57,:);
v17 = X(58,:);
v18 = X(59,:);
v19 = X(60,:);
v20 = X(61,:);
v21 = X(62,:);
v22 = X(63,:);
v23 = X(64,:);
v24 = X(65,:);
v25 = X(66,:);
v26 = X(67,:);
v27 = X(68,:);
v28 = X(69,:);
v29 = X(70,:);
v30 = X(71,:);



% Netzwerk
p1 = X(72,:);
p2 = X(73,:);
p3 = X(74,:);
p4 = X(75,:);
p5 = X(76,:);
p6 = X(77,:);
p7 = X(78,:);
p8 = X(79,:);
p9 = X(80,:);
p10 = X(81,:);
p11 = X(82,:);
p12 = X(83,:);
p13 = X(84,:);
p14 = X(85,:);
p15 = X(86,:);
p16 = X(87,:);
p17 = X(88,:);
p18 = X(89,:);
p19 = X(90,:);
p20 = X(91,:);
p21 = X(92,:);
p22 = X(93,:);
p23 = X(94,:);
p24 = X(95,:);
p25 = X(96,:);
p26 = X(97,:);
p27 = X(98,:);
p28 = X(99,:);
p29 = X(100,:);
p30 = X(101,:);


q1 = X(102,:);
q2 = X(103,:);
q3 = X(104,:);
q4 = X(105,:);
q5 = X(106,:);
q6 = X(107,:);
q7 = X(108,:);
q8 = X(109,:);
q9 = X(110,:);
q10 = X(111,:);
q11 = X(112,:);
q12 = X(113,:);
q13 = X(114,:);
q14 = X(115,:);
q15 = X(116,:);
q16 = X(117,:);
q17 = X(118,:);
q18 = X(119,:);
q19 = X(120,:);
q20 = X(121,:);
q21 = X(122,:);
q22 = X(123,:);
q23 = X(124,:);
q24 = X(125,:);
q25 = X(126,:);
q26 = X(127,:);
q27 = X(128,:);
q28 = X(129,:);
q29 = X(130,:);
q30 = X(131,:);


%v9 = X(42,:);


%% Vektoren f??r 9 Bus Netzwerk
V =  [ V1 V2 v3 v4 V5 v6 v7 V8 v9 v10 V11 v12 V13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 ].';

Pw = [ p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30].';
Q =  [ q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16 q17 q18 q19 q20 q21 q22 q23 q24 q25 q26 q27 q28 q29 q30].';

Theta = [0 theata2 theta3 theta4 theata5 theta6 theta7 theata8 theta9 theta10 theata11 theta12 theata13 theta14 theta15 ...
           theta16 theta17 theta18 theta19 theta20 theta21 theta22 theta23 theta24 theta25 theta26 theta27 theta28 theta29 theta30];

n = length (Pw);
Theta_hk = Theta'*ones(1,n)-ones(n,1)*Theta ;

Pe1  = (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-I.theta1));
Pe2  = (I.E2*V2*abs(P.xd2)*cos(angle(P.xd2)+delta2-theata2));
Pe5  = (I.E5*V5*abs(P.xd5)*cos(angle(P.xd5)+delta5-theata5));
Pe8  = (I.E8*V8*abs(P.xd8)*cos(angle(P.xd8)+delta8-theata8));
Pe11 = (I.E11*V11*abs(P.xd11)*cos(angle(P.xd11)+delta11-theata11));
Pe13 = (I.E13*V13*abs(P.xd13)*cos(angle(P.xd13)+delta13-theata13));

out = [
P.omegaS*omega1
(1/P.M1)*(I.Pm1-Pe1 - (P.D1*omega1))

P.omegaS*omega2
(1/P.M2)*(I.Pm2-Pe2 - (P.D2*omega2))

P.omegaS*omega5
(1/P.M5)*(I.Pm5-Pe5 - (P.D5*omega5))

P.omegaS*omega8
(1/P.M8)*(I.Pm8-Pe8 - (P.D8*omega8))

P.omegaS*omega11
(1/P.M11)*(I.Pm11-Pe11 - (P.D11*omega11))

P.omegaS*omega13
(1/P.M13)*(I.Pm13-Pe13 - (P.D13*omega13))

p1 - (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1));
p2 - (I.E2*V2*abs(P.xd2)*cos(angle(P.xd2)+delta2-theata2));
p3 - I.p3 * v3^2 / I.v3^2;
p4 - I.p4 * v4^2 / I.v4^2;
p5 - (I.E5*V5*abs(P.xd5)*cos(angle(P.xd5)+delta5-theata5));
p6 - I.p6 * v6^2 / I.v6^2;
p7 - I.p7 * v7^2 / I.v7^2;
p8 - (I.E8*V8*abs(P.xd8)*cos(angle(P.xd8)+delta8-theata8));
p9 - I.p9 * v9^2 / I.v9^2;
p10 - I.p10 * v10^2 / I.v10^2;
p11 - (I.E11*V11*abs(P.xd11)*cos(angle(P.xd11)+delta11-theata11));
p12 - I.p12 * v12^2 / I.v12^2;
p13 - (I.E13*V13*abs(P.xd13)*cos(angle(P.xd13)+delta13-theata13));
p14 - I.p14 * v14^2 / I.v14^2;
p15 - I.p15 * v15^2 / I.v15^2;
p16 - I.p16 * v16^2 / I.v16^2;
p17 - I.p17 * v17^2 / I.v17^2;
p18 - I.p18 * v18^2 / I.v18^2;
p19 - I.p19 * v19^2 / I.v19^2;
p20 - I.p20 * v20^2 / I.v20^2;
p21 - I.p21 * v21^2 / I.v21^2;
p22 - I.p22 * v22^2 / I.v22^2;
p23 - I.p23 * v23^2 / I.v23^2;
p24 - I.p24 * v24^2 / I.v24^2;
p25 - I.p25 * v25^2 / I.v25^2;
p26 - I.p26 * v26^2 / I.v26^2;
p27 - I.p27 * v27^2 / I.v27^2;
p28 - I.p28 * v28^2 / I.v28^2;
p29 - I.p29 * v29^2 / I.v29^2;
p30 - I.p30 * v30^2 / I.v30^2;


%q1 + (I.E1*V1*abs(P.xd1)*sin(angle(P.xd1)+delta1)) - (V1^2*abs(P.xd1)*sin(angle(P.xd1)));
q2 + (I.E2*V2*abs(P.xd2)*sin(angle(P.xd2)+delta2-theata2)) - (V2^2*abs(P.xd2)*sin(angle(P.xd2)));
q3 - I.q3 * v3^2 / I.v3^2;
q4 - I.q4 * v4^2 / I.v4^2;
q5 + (I.E5*V5*abs(P.xd5)*sin(angle(P.xd5)+delta5-theata5)) - (V5^2*abs(P.xd5)*sin(angle(P.xd5)));
q6 - I.q6 * v6^2 / I.v6^2;
q7 - I.q7 * v7^2 / I.v7^2;
q8 + (I.E8*V8*abs(P.xd8)*sin(angle(P.xd8)+delta8-theata8)) - (V8^2*abs(P.xd8)*sin(angle(P.xd8)));
q9 - I.q9 * v9^2 / I.v9^2;
q10 - I.q10 * v10^2 / I.v10^2;
q11 + (I.E11*V11*abs(P.xd11)*sin(angle(P.xd11)+delta11-theata11)) - (V11^2*abs(P.xd11)*sin(angle(P.xd11)));
q12 - I.q12 * v12^2 / I.v12^2;
q13 + (I.E13*V13*abs(P.xd13)*sin(angle(P.xd13)+delta13-theata13)) - (V13^2*abs(P.xd13)*sin(angle(P.xd13)));
q14 - I.q14 * v14^2 / I.v14^2;
q15 - I.q15 * v15^2 / I.v15^2;
q16 - I.q16 * v16^2 / I.v16^2;
q17 - I.q17 * v17^2 / I.v17^2;
q18 - I.q18 * v18^2 / I.v18^2;
q19 - I.q19 * v19^2 / I.v19^2;
q20 - I.q20 * v20^2 / I.v20^2;
q21 - I.q21 * v21^2 / I.v21^2;
q22 - I.q22 * v22^2 / I.v22^2;
q23 - I.q23 * v23^2 / I.v23^2;
q24 - I.q24 * v24^2 / I.v24^2;
q25 - I.q25 * v25^2 / I.v25^2;
q26 - I.q26 * v26^2 / I.v26^2;
q27 - I.q27 * v27^2 / I.v27^2;
q28 - I.q28 * v28^2 / I.v28^2;
q29 - I.q29 * v29^2 / I.v29^2;
q30 - I.q30 * v30^2 / I.v30^2;

V.*( (real(xs).*cos(Theta_hk) + imag(xs).*sin(Theta_hk))* V ) - Pw; 
V.*( (real(xs).*sin(Theta_hk) - imag(xs).*cos(Theta_hk))* V ) - Q; 
];
