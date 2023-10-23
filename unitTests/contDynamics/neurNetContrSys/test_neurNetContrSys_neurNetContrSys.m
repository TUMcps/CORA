function res = test_neurNetContrSys_neurNetContrSys
% test_neurNetContrSys_neurNetContrSys - unit test function of 
%    neurNetContrSys constructor
%
%
% Syntax:
%    res = test_neurNetContrSys_neurNetContrSys
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Tobias Ladner
% Written:       23-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% contDynamics
f = @(x,u) [x(2) + u(2); (1-x(1)^2)*x(2) - x(1) + u(1)];
sysOL = nonlinearSys(f);

% neural network controller
W1 = rand(100,2); b1 = rand(100,1);
W2 = rand(1,100); b2 = rand(1,1);
nn = neuralNetwork({ ...
    nnLinearLayer(W1, b1); ...
    nnReLULayer(); ...
    nnLinearLayer(W2, b2); ...
    nnReLULayer(); ...
});

% neural network controlled system
dt = 0.01;
sys = neurNetContrSys(sysOL,nn,dt);

% check old neural network input
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

% check wrong input
try    
    sys = neurNetContrSys();
    % should have thrown an error
    res = false;
end
try    
    sys = neurNetContrSys(sysOL,nn,[0.1, 0.2]);
    % should have thrown an error
    res = false;
catch
end

try
    W1 = rand(100, 3);
    nn = neuralNetwork({ ...
        nnLinearLayer(W1, b1); ...
        nnReLULayer(); ...
        nnLinearLayer(W2, b2); ...
        nnReLULayer(); ...
    });
    sys = neurNetContrSys(sysOL,nn,[0.1, 0.2]);
    res = false;
catch
end

res = true;

% ------------------------------ END OF CODE ------------------------------
