function res = test_nnActivationLayer_instantiateFromString()
% test_nnActivationLayer_instantiateFromString - tests the instantiation
%    of a activation layer from a string
%
% Syntax:  
%    res = test_nnActivationLayer_instantiateFromString()
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

% Author:       Tobias Ladner
% Written:      28-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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

try
    % should throw an error
    layer = nnActivationLayer.instantiateFromString('unkown');
    res = false;
end


%------------- END OF CODE --------------
