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
% Last update:   20-February-2026 (TL, add linearSys/linearSysDT/nonlinearSysDT/linParamSys)
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

% linearSys ---

A = [-1 1; 0 -2]; B = [0; 1];
sysOL = linearSys(A, B);

W1 = rand(10,2); b1 = rand(10,1);
W2 = rand(1,10); b2 = rand(1,1);
nn = neuralNetwork({ ...
    nnLinearLayer(W1, b1); nnReLULayer(); ...
    nnLinearLayer(W2, b2); nnReLULayer(); ...
});

dt = 0.05;
sys = neurNetContrSys(sysOL, nn, dt);

% linearSysDT ---

A = [0.9 0.1; 0 0.8]; B = [0; 1];
sysOL = linearSysDT(A, B, dt);

sys = neurNetContrSys(sysOL, nn, dt);

% nonlinearSysDT ---

f = @(x,u) [x(1) + 0.05*x(2); x(2) + 0.05*u(1)];
sysOL = nonlinearSysDT(f, dt);

sys = neurNetContrSys(sysOL, nn, dt);

% linParamSys (numeric A) ---

Ac = [-1 1; 0 -2]; Aw = [0.1 0; 0 0.1];
A = intervalMatrix(Ac, Aw); B = [0; 1];
sysOL = linParamSys(A, B, 'varParam');

sys = neurNetContrSys(sysOL, nn, dt);

% check wrong input
assertThrowsAs(@neurNetContrSys,'MATLAB:structRefFromNonStruct');
assertThrowsAs(@neurNetContrSys,'CORA:wrongValue',sysOL,nn,[0.1, 0.2]);

W1 = rand(100,3);
b1 = rand(100,1);
nn = neuralNetwork({ ...
    nnLinearLayer(W1, b1); ...
    nnReLULayer(); ...
    nnLinearLayer(W2, b2); ...
    nnReLULayer(); ...
});
assertThrowsAs(@neurNetContrSys,'CORA:wrongValue',sysOL,nn,[0.1, 0.2]);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
