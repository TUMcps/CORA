%%
%parameter update
% func = @(x,u)vehicle(x,u);
sys = nonlinearSys(@vehicle, 11, 1);
% parameter
params.tFinal = 0.05;
% params.R0 = zonotope([0.015; 0.015; 0.015; 0.015; 0.015; 0.015; 0.015; 0.015; 0.015; 0.015; 0.015],diag([0.005; 0.005; 0.005; 0.005; 0.005; 0.005; 0.005; 0.005; 0.005; 0.005; 0.005]));
R0 = interval([0.1; 0.1; 0; pi/6; pi/6; pi/6; pi/6; 0.1; 0.1; 0.1; 0.1], [0.15; 0.15; 0; pi/5; pi/5; pi/5; pi/5; 0.15; 0.15; 0.15; 0.15]);
params.R0 = zonotope(R0);
% reachability settings
% options.algInner = 'parallelo';
options.algInner = 'scale';
options.alg = 'lin-adaptive';
options.verbose = true;
% options.alg = 'poly';

% options.timeStep = params.tFinal / 500;
% options.taylorTerms = 10;
% options.zonotopeOrder = 100;
% options.alg = 'lin';
% options.tensorOrder = 2;
% options.errorOrder = 10;
% options.intermediateOrder = 20;

% reachability analysis
R = reach(sys,params,options);

% Rin = reachInner(sys,params,options);