function res = test_nn_nnLipConstrLinearLayer()
% test_nn_nnLipConstrLinearLayer - tests constructor of nnLipConstrLinearLayer
%
% Syntax:
%    res = test_nn_nnLipConstrLinearLayer()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Lukas Koller
% Written:       18-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Instantiate a linear layer with fixed weights.
W = rand(5,3); 
b = rand(5,1);
lambda = 3;
name = "TestLayer";
layer = nnLipConstrLinearLayer(W, b, lambda, name);

% Check if the attributes are correctly assigned.
assert(compareMatrices(W, layer.W))
assert(compareMatrices(b, layer.b))
assert(lambda == layer.lambda)
assert(strcmp(name,layer.name))

% Check that the normalization.
options.nn.train.backprop = true;
layer.normWeights(options);
% Obtain the normalization matrix.
W_norm = layer.backprop.store.W_norm;
assert(compareMatrices(W*W_norm,layer.W));
% Normalizing twice should not do anything.
layer.normWeights(options);
assert(compareMatrices(W*W_norm,layer.W));
assert(compareMatrices(layer.backprop.store.W_norm,eye(3)));

% Check variable input.
layer = nnLipConstrLinearLayer(W);
assert(sum(layer.b) == 0)

% Check wrong input.
assertThrowsAs(@nnLipConstrLinearLayer,'MATLAB:minrhs');

% Check for a dimension missmatch.
W = rand(4,3); 
b = rand(10,1);
assertThrowsAs(@nnLipConstrLinearLayer,'CORA:wrongInputInConstructor',W,b);

% Test completed.
res = true;

% ------------------------------ END OF CODE ------------------------------
