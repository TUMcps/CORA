function layer = instantiateFromString(activation)
% instantiateFromString - creates an activation layer from string
%
% Syntax:
%    layer = nnActivationLayer.instantiateFromString(activation)
%
% Inputs:
%    activation - string, one of {'relu','sigmoid','tanh', 'softmax', 
%       'identity', 'none'}
%
% Outputs:
%    layer - nnActivationLayer
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnLayer

% Authors:       Tobias Ladner
% Written:       24-June-2022
% Last update:   02-November-2022 (lower)
%                29-November-2022 (switch, softmax)
%                30-November-2022 (none, identity)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
activation = lower(activation);
possibleActivations = {'relu','sigmoid','tanh','softmax', 'identity', 'none'};
inputArgsCheck({{activation,'str', possibleActivations}});

switch activation
    case "relu"
        layer = nnReLULayer();
    case "sigmoid"
        layer = nnSigmoidLayer();
    case "tanh"
        layer = nnTanhLayer();
    case "softmax"
        layer = nnSoftmaxLayer();
    case "identity"
        layer = nnIdentityLayer();
    case "none"
        layer = nnIdentityLayer();
    otherwise
        % should not be executed anyway due to inputArgsCheck. 
        throw(CORAerror('CORA:wrongValue', 'first', ...
            strjoin(possibleActivations, ', ')));
end

end

% ------------------------------ END OF CODE ------------------------------
