function text = benchmark_hybrid_reach_ARCH23_powerTrain_DTN06()
% benchmark_hybrid_reach_ARCH23_powerTrain_DTN06 - powertrain example 
%    (see [1]) from the 2023 ARCH competition (Test case DTN06)
% 
%
% Syntax:
%    benchmark_hybrid_reach_ARCH23_powerTrain_DTN06()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] M. Althoff and B. H. Krogh. "Avoiding geometric intersection 
%       operations in reachability analysis of hybrid systems." HSCC 2012.

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       21-September-2011
% Last update:   12-March-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% System Dynamics ---------------------------------------------------------

% get system dynmics
[HA, Zcenter, Zdelta, B1, B2, B3, c1, c2, c3] = initPowerTrain(51); 

% split initial set
R01 = zonotope(Zcenter + 0.75*Zdelta,0.25*Zdelta);
R02 = zonotope(Zcenter + 0.25*Zdelta,0.25*Zdelta);
R03 = zonotope(Zcenter - 0.25*Zdelta,0.25*Zdelta);
R04 = zonotope(Zcenter - 0.75*Zdelta,0.25*Zdelta);


% Parameters --------------------------------------------------------------

% initial set, initial location, and final time
params.startLoc = 3;
params.tFinal = 0.2; 

% system input
u = [-5; 0];            % acceleration and load torque

params.U{1} = zonotope(B1*u + c1);
params.U{2} = zonotope(B2*u + c2); 
params.U{3} = zonotope(B3*u + c3);


% Reachability Settings ---------------------------------------------------

% continuous reachability
options.taylorTerms = 20;
options.zonotopeOrder = 20;

options.timeStep{1} = 2e-4;
options.timeStep{2} = 1e-4;
options.timeStep{3} = 1e-4;

% intersection with guard sets
options.guardIntersect = 'zonoGirard';
options.enclose = {'pca'};


% Simulation 1 ------------------------------------------------------------

% simulation parameter
paramsSim = params;
paramsSim = rmfield(paramsSim,'U');

paramsSim.u{1} = center(params.U{1});
paramsSim.u{2} = center(params.U{2});
paramsSim.u{3} = center(params.U{3});

% obtain random simulation results
N = 10;
simPoints = cell(N,1);
simRes = [];

for i = 1:N
    
    % set initial state
    if i == 1
        paramsSim.x0 = Zcenter + Zdelta;
    elseif i == 2
        paramsSim.x0 = Zcenter - Zdelta;
    else
        paramsSim.x0 = randPoint(zonotope(Zcenter,Zdelta));
    end 
    
    % simulate hybrid automaton
    [t,x,loc] = simulate(HA,paramsSim); 
    
    % store results
    simPoints{i} = x{end}(end,:)';
    simRes = add(simRes,simResult(x,t,loc)); 
end


% Reachability Analysis 1 -------------------------------------------------

tStart = tic;

params.R0 = R01;
R_11 = reach(HA,params,options);

params.R0 = R02;
R_12 = reach(HA,params,options);

params.R0 = R03;
R_13 = reach(HA,params,options);

params.R0 = R04;
R_14 = reach(HA,params,options);

tComp1 = toc(tStart);

R = add(add(add(R_11,R_12),R_13),R_14);


% Parameter ---------------------------------------------------------------

% update initial set
R01 = R_11.timePoint.set{end};
R02 = R_12.timePoint.set{end};
R03 = R_13.timePoint.set{end};
R04 = R_14.timePoint.set{end};

% update start and final time
params.tStart = 0.2;
params.tFinal = 2;

% update system input
u = [5; 0];             % acceleration and load torque

params.U{1} = zonotope(B1*u + c1);
params.U{2} = zonotope(B2*u + c2); 
params.U{3} = zonotope(B3*u + c3);


% Simulation 2 ------------------------------------------------------------

% simulation parameter
paramsSim = params;
paramsSim = rmfield(paramsSim,'U');
paramsSim = rmfield(paramsSim,'R0');

paramsSim.u{1} = center(params.U{1});
paramsSim.u{2} = center(params.U{2});
paramsSim.u{3} = center(params.U{3});

% obtain random simulation results
for i = 1:N
    
    % set initial state
    paramsSim.x0 = simPoints{i};
    
    % simulate hybrid automaton
    [t,x,loc] = simulate(HA,paramsSim); 
    
    % store resulets
    simRes = add(simRes,simResult(x,t,loc));
end


% Reachability Analysis 2 -------------------------------------------------

tStart = tic;

params.R0 = R01;
R_21 = reach(HA,params,options);

params.R0 = R02;
R_22 = reach(HA,params,options);

params.R0 = R03;
R_23 = reach(HA,params,options);

params.R0 = R04;
R_24 = reach(HA,params,options);

tComp2 = toc(tStart);

disp(['computation time: ',num2str(tComp1 + tComp2)]);

R = add(add(add(add(R_21,R_22),R_23),R_24),R);


% Visualization -----------------------------------------------------------

figure; hold on; box on;
useCORAcolors('CORA:contDynamics')

% plot reachable set
plot(R,[1,3]);
plot(R(1).R0,[1,3]);

% plot simulation
plot(simRes,[1,3]);

% guard set
plot(specification(interval([0.03; -10], [0.03; 90])));
plot(specification(interval([-0.03; -10], [-0.03; 90])));

% formatting
xlabel('x_1');
ylabel('x_3');
xlim([-0.1,0.2]);
ylim([-10,90]);

% example completed
res = true;
text = ['Powertrain,DTN06,',num2str(res),',',num2str(tComp1+tComp2)];

% ------------------------------ END OF CODE ------------------------------
