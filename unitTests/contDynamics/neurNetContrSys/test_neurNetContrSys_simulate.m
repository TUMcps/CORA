function res = test_neurNetContrSys_simulate
% test_neurNetContrSys_simulate - unit test function of 
%    simulate function of neurNetContrSys
%
%
% Syntax:
%    res = test_neurNetContrSys_simulate
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Tobias Ladner
% Written:       28-November-2022
% Last update:   20-February-2026 (TL, add linearSys/linearSysDT/nonlinearSysDT/linParamSys)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dynamic system
f = @(x,u) [x(2) + u(2); (1-x(1)^2)*x(2) - x(1) + u(1)];
sysOL = nonlinearSys(f);

% neural network controller
layers = cell(4, 1);
W1 = rand(10,2); b1 = rand(10,1);
layers{1} = nnLinearLayer(W1, b1);
layers{2} = nnSigmoidLayer();
W2 = rand(2,10); b2 = rand(2,1);
layers{3} = nnLinearLayer(W2, b2);
layers{4} = nnSigmoidLayer();
nn = neuralNetwork(layers);

% neural network controlled system
dt = 0.01;
sys = neurNetContrSys(sysOL,nn,dt);

params.x0 = [1;2];
params.tFinal = 1;

[t,x] = simulate(sys,params);

figure
plot(x(1,:),x(2,:),'k', 'DisplayName', 'Simulation');
legend()
close(gcf);

% shared nn (2 inputs -> 1 output) and short params for remaining systems
W1 = rand(10,2); b1 = rand(10,1);
W2 = rand(1,10); b2 = rand(1,1);
nn = neuralNetwork({ ...
    nnLinearLayer(W1, b1); nnReLULayer(); ...
    nnLinearLayer(W2, b2); nnReLULayer(); ...
});
dt = 0.05;
params_short.x0 = [1; 0];
params_short.tFinal = 0.2;

% linearSys ---

A = [-1 1; 0 -2]; B = [0; 1];
sysOL = linearSys(A, B);
sys = neurNetContrSys(sysOL, nn, dt);
simulate(sys, params_short);

% linearSysDT ---

A = [0.9 0.1; 0 0.8]; B = [0; 1];
sysOL = linearSysDT(A, B, dt);
sys = neurNetContrSys(sysOL, nn, dt);
simulate(sys, params_short);

% nonlinearSysDT ---

f = @(x,u) [x(1) + 0.05*x(2); x(2) + 0.05*u(1)];
sysOL = nonlinearSysDT(f, dt);
sys = neurNetContrSys(sysOL, nn, dt);
simulate(sys, params_short);

% linParamSys ---

Ac = [-1 1; 0 -2]; Aw = [0.1 0; 0 0.1];
A = intervalMatrix(Ac, Aw); B = [0; 1];
sysOL = linParamSys(A, B, 'varParam');
sys = neurNetContrSys(sysOL, nn, dt);
simulate(sys, params_short);

res = true;

% ------------------------------ END OF CODE ------------------------------
