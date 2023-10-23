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
% Last update:   ---
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
plot(x(:,1),x(:,2),'k', 'DisplayName', 'Simulation');
legend()
close(gcf);

res = true;

% ------------------------------ END OF CODE ------------------------------
