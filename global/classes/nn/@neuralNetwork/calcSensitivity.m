function S = calcSensitivity(obj, x, varargin)
% calcSensitivity - calculates input-output sensitivity matrix at x
%       rows correspond to output neurons, columns to input neurons
%       sensitivity of layer i will be stored in obj.layers{i}.sensitivity
%
% Syntax:
%    S = calcSensitivity(obj, x)
%
% Inputs:
%    obj - object of class neuralNetwork
%    x - point from input space
%    evParams
%
% Outputs:
%    S - sensitivity matrix at x
%
% References:
%    [1] Zurada, J. M., et al. "Sensitivity analysis for minimization of
%           input data dimension for feedforward neural network"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork, nnHelper/validateEvaluateParams

% Authors:       Tobias Ladner
% Written:       14-April-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 2
    throw(CORAerror('CORA:notEnoughInputArgs', 2))
elseif nargin > 3
    throw(CORAerror("CORA:tooManyInputArgs", 3))
end
evParams = setDefaultValues({struct}, varargin);
evParams = nnHelper.validateEvaluateParams(evParams);

% forward propagation
xs = cell(length(obj.layers), 1);
for i = 1:length(obj.layers)
    xs{i} = x;
    layer_i = obj.layers{i};
    x = layer_i.evaluateNumeric(x, evParams);
end

% calculate sensitivity [1]
S = 1;

% backward propagation
for i = length(obj.layers):-1:1
    layer_i = obj.layers{i};
    S = layer_i.evaluateSensitivity(S, xs{i}, evParams);

    % save sensitivity at layer i for refinement
    layer_i.sensitivity = S;
end

end

% ------------------------------ END OF CODE ------------------------------
