function text = benchmark_hybrid_reach_ARCH23_powerTrain_DTN05()
% benchmark_hybrid_reach_ARCH23_powerTrain_DTN05 - powertrain example 
%    (see [1]) from the 2023 ARCH competition (Test case DTN05)
% 
% Syntax:
%    benchmark_hybrid_reach_ARCH23_powerTrain_DTN05()
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

[HA, Zcenter, Zdelta, B1, B2, B3, c1, c2, c3] = initPowerTrain(51); 
Zdelta = 0.3*Zdelta;


% Parameters --------------------------------------------------------------

% initial set, initial location, and final time
params.R0 = zonotope([Zcenter,Zdelta]);
params.startLoc = 3;
params.tFinal = 0.2; 

% system input
u = [-5; 0];            % acceleration and load torque

params.U{1} = B1*zonotope(0) + B1*u + c1;
params.U{2} = B2*zonotope(0) + B2*u + c2; 
params.U{3} = B3*zonotope(0) + B3*u + c3;


% Reachability Settings ---------------------------------------------------

% continuous reachability
options.taylorTerms = 20;
options.zonotopeOrder = 20;

options.timeStep = 1e-4;

% intersection with guard sets
options.guardIntersect = 'zonoGirard';
options.enclose = {'pca'};


% Simulation 1 ------------------------------------------------------------

% simulation parameter
paramsSim = params;
paramsSim = rmfield(paramsSim,'U');
paramsSim = rmfield(paramsSim,'R0');

paramsSim.u{1} = center(params.U{1});
paramsSim.u{2} = center(params.U{2});
paramsSim.u{3} = center(params.U{3});

paramsSim.startLoc = params.startLoc; 

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
        paramsSim.x0 = randPoint(params.R0);
    end 
    
    % simulate hybrid automaton
    [t,x,loc] = simulate(HA,paramsSim); 
    
    % store results
    simPoints{i} = x{end}(end,:)';
    simRes = add(simRes,simResult(x,t,loc)); 
end


% Reachability Analysis 1 -------------------------------------------------

tStart = tic;
[R,res1] = reach(HA,params,options);
tComp1 = toc(tStart);


% Parameter ---------------------------------------------------------------

% update initial set
params.R0 = R.timePoint.set{end};

% update start and final time
params.tStart = 0.2;
params.tFinal = 2;

% update system input
u = [5; 0];             % acceleration and load torque

params.U{1} = B1*zonotope(0) + B1*u + c1;
params.U{2} = B2*zonotope(0) + B2*u + c2; 
params.U{3} = B3*zonotope(0) + B3*u + c3;


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
[Rtemp,res2] = reach(HA,params,options);
tComp2 = toc(tStart);

R = add(R,Rtemp);

disp(['computation time: ',num2str(tComp1 + tComp2)]);


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
text = ['Powertrain,DTN05,',num2str(res1&res2),',',num2str(tComp1+tComp2)];

% ------------------------------ END OF CODE ------------------------------
