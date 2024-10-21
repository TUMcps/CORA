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

activation = 'ReLU';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnActivationLayer'));

activation = 'tanh';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnTanhLayer'));

activation = 'sigmoid';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnSigmoidLayer'));

activation = 'softmax';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnSoftmaxLayer'));

activation = 'identity';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnIdentityLayer'));

activation = 'none';
layer = nnActivationLayer.instantiateFromString(activation);
assert(isa(layer, 'nnIdentityLayer'));

assertThrowsAs(@nnActivationLayer.instantiateFromString,'CORA:wrongValue','unknown');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
