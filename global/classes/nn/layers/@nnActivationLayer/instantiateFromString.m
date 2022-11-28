function layer = instantiateFromString(activation)
% instantiateFromString - creates an activation layer from string
%
% Syntax:
%    layer = nnActivationLayer.instantiateFromString(activation)
%
% Inputs:
%    activation - one of {'ReLU','sigmoid','tanh'}
%
% Outputs:
%    layer - nnActivationLayer
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnLayer

% Author:        Tobias Ladner
% Written:       24-June-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% check input arguments
inputArgsCheck({{activation,'str',{'ReLU','sigmoid','tanh'}}});

if strcmp(activation, "ReLU")
    layer = nnReLULayer();
elseif strcmp(activation, "sigmoid")
    layer = nnSigmoidLayer();
elseif strcmp(activation, "tanh")
    layer = nnTanhLayer();
end

%------------- END OF CODE --------------