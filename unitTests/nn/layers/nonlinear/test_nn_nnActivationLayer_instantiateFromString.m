function res = test_nn_nnActivationLayer_instantiateFromString()
% test_nn_nnActivationLayer_instantiateFromString - tests the instantiation
%    of a activation layer from a string
%
% Syntax:
%    res = test_nn_nnActivationLayer_instantiateFromString()
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

% Authors:       Tobias Ladner
% Written:       28-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% relu
activation = 'ReLU';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnActivationLayer'));

% tanh
activation = 'tanh';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnTanhLayer'));

% sigmoid
activation = 'sigmoid';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnSigmoidLayer'));

% softmax
activation = 'softmax';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnSoftmaxLayer'));

% identity
activation = 'identity';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnIdentityLayer'));

% none
activation = 'none';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnIdentityLayer'));

% no activation given
assertThrowsAs(@nnActivationLayer.instantiateFromString,'CORA:wrongValue','unknown');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
