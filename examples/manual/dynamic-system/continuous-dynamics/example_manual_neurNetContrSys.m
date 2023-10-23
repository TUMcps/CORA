function example_manual_neurNetContrSys()
% example_manual_neurNetContrSys - example from the manual demonstrating the 
% neurNetContrSys constructor as defined in the manual
%
% Syntax:
%   example_manual_neurNetContrSys()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:        Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% open-loop system
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
})

% neural network controlled system
dt = 0.01;
sys = neurNetContrSys(sysOL,nn,dt);

% ------------------------------ END OF CODE ------------------------------
