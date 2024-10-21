function res = example_manual_neuralNetwork()
% example_manual_neuralNetwork - example from the manual demontrating the
% neuralNetwork class as defined in the manual
%
% Syntax:
%   example_manual_neuralNetwork()
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

% Authors:       Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init weight and bias
W1 = rand(3,2); b1 = rand(3,1);
W2 = rand(2,3); b2 = rand(2,1);

% neural network
nn = neuralNetwork({ ...
    nnLinearLayer(W1, b1); ...
    nnReLULayer(); ...
    nnLinearLayer(W2, b2); ...
    nnReLULayer(); ...
});

 % input set
 c = [4;4];
 G = [2 1 2; 0 2 2];
 E = [1 0 3;0 1 1];
 GI = [];
 X = polyZonotope(c,G,GI,E);
 
 % settings
 options = struct;
 options.nn.poly_method = 'regression';
 options.nn.num_generators = 1000;
 
 % evaluation
 Y = nn.evaluate(X, options);

 % example completed
 res = true;

% ------------------------------ END OF CODE ------------------------------
