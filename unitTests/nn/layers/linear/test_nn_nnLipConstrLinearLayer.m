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

% Specify number of input and output neurons.
nIn = 31;
nOut = 27;

% Instantiate a linear layer with fixed weights.
W = rand([nOut nIn]); 
b = zeros([nOut 1]); % Bias does not matter.
lambda = randi([1 5]);
name = 'TestLayer';
layer = nnLipConstrLinearLayer(W, b, lambda, name);

% Check if the attributes are correctly assigned.
assert(compareMatrices(W, layer.W))
assert(compareMatrices(b, layer.b))
assert(lambda == layer.lambda)
assert(strcmp(name,layer.name))

% Check that the normalization.
% Specify number of inputs.
N = 10;
% Generate a random inputs.
x = rand([nIn N]);
% Compute the output.
y = layer.evaluate(x);
% Compute a lower bound of the Lipschitz constant from the input and
% outputs.
dy = pdist(y','chebychev'); % Compute l-inf norm between outputs
dx = pdist(x','cityblock'); % COmpute l-1 norm between inputs
% Compute all possible lowerbounds for Lipschitz constants.
labmdas = dy./dx;
% Check lower bounds.
assert(all(lambda >= labmdas | isnan(lambda),'all'));

% Check bias.
b = rand([nOut 1]);
layer = nnLipConstrLinearLayer(zeros([nOut nIn]),b);
% Compute the outputs.
y = layer.evaluate(x);
% Check the outputs.
assert(all(y == b,'all'));

% Check variable input.
layer = nnLipConstrLinearLayer(W);
assert(sum(layer.b) == 0)

% Check wrong input.
assertThrowsAs(@nnLipConstrLinearLayer,'MATLAB:minrhs');

% Check for a dimension missmatch.
W = rand([nOut nIn]); 
b = rand([nOut+3 1]);
assertThrowsAs(@nnLipConstrLinearLayer,'CORA:wrongInputInConstructor',W,b);

% Test completed.
res = true;

% ------------------------------ END OF CODE ------------------------------
