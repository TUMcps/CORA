function res = example_simulate()

% systems
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;
samplingTime = 0.02;

% continuous
fiveDimSys = linearSys('fiveDimSys',A,B);
% discrete
fiveDimSys_DT = linearSysDT(fiveDimSys,samplingTime);

% parameters
params.tFinal = 5;
params.R0 = zonotope([ones(5,1),0.1*diag(ones(5,1))]);
params.U = zonotope(interval([0.9; -0.25; -0.1; 0.25; -0.75], ...
                             [1.1; 0.25; 0.1; 0.75; -0.25]));

% zonotope order
options.zonotopeOrder = 20; 

% reachability analysis
R_DT = reach(fiveDimSys_DT, params, options);

% additional options for continuous systems
options.timeStep = samplingTime;
options.taylorTerms = 4;
R = reach(fiveDimSys, params, options);

% simulate
params.x0 = randPoint(params.R0);
params.u = randPoint(params.U);
[t,x] = simulate(fiveDimSys,params);
[t,x_DT] = simulate(fiveDimSys_DT,params);

% plot
figure;
subplot(1,2,1); hold on;
title("Continuous");
hR = plot(R);
hSim = scatter(x(:,1),x(:,2),8,'r','filled');
legend([hR,hSim],'reach','simulate');

subplot(1,2,2); hold on;
title("Discrete");
hR = plot(R_DT);
hSim = scatter(x_DT(:,1),x_DT(:,2),8,'r','filled');
legend([hR,hSim],'reach','simulate');

% example completed
res = 1;

%------------- END OF CODE --------------