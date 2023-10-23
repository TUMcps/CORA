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
% See also: -

% Authors:       Tobias Ladner
% Written:       28-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

activation = 'ReLU';
layer = nnActivationLayer.instantiateFromString(activation);
res = res & isa(layer, 'nnActivationLayer');

activation = 'tanh';
layer = nnActivationLayer.instantiateFromString(activation);
res = res & isa(layer, 'nnTanhLayer');

activation = 'sigmoid';
layer = nnActivationLayer.instantiateFromString(activation);
res = res & isa(layer, 'nnSigmoidLayer');

activation = 'softmax';
layer = nnActivationLayer.instantiateFromString(activation);
res = res & isa(layer, 'nnSoftmaxLayer');

activation = 'identity';
layer = nnActivationLayer.instantiateFromString(activation);
res = res & isa(layer, 'nnIdentityLayer');

activation = 'none';
layer = nnActivationLayer.instantiateFromString(activation);
res = res & isa(layer, 'nnIdentityLayer');

try
    % should throw an error
    layer = nnActivationLayer.instantiateFromString('unknown');
    res = false;
end


% ------------------------------ END OF CODE ------------------------------
