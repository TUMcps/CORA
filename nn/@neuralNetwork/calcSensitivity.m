function [S,y] = calcSensitivity(obj, x, varargin)
% calcSensitivity - calculates input-output sensitivity matrix at x
%       rows correspond to output neurons, columns to input neurons
%       sensitivity of layer i will be stored in obj.layers{i}.sensitivity
%
% Syntax:
%    [S,y] = calcSensitivity(obj, x)
%
% Inputs:
%    obj - object of class neuralNetwork
%    x - point from input space
%    options - options for neural network evaluation (stored in options.nn)
%    store_sensitivty - {0,1} if sensitivity should be stored in each
%       layer; default: 1
%
% Outputs:
%    S - sensitivity matrix at x
%    y - output of the neural network for x
%
% References:
%    [1] Zurada, J. M., et al. "Sensitivity analysis for minimization of
%           input data dimension for feedforward neural network"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork, nnHelper/validateNNoptions

% Authors:       Tobias Ladner, Lukas Koller
% Written:       14-April-2022
% Last update:   16-January-2024 (return neural network output)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(2,4);
[options,storeSensitivty] = setDefaultValues({struct,true}, varargin);
options = nnHelper.validateNNoptions(options);

% forward propagation
xs = cell(length(obj.layers), 1);
for i = 1:length(obj.layers)
    xs{i} = x;
    layer_i = obj.layers{i};
    x = layer_i.evaluateNumeric(x, options);
end
y = x;

% calculate sensitivity [1]
S = ones([1 1 size(y,2)],'like',y);

% backward propagation
for i = length(obj.layers):-1:1
    layer_i = obj.layers{i};
    S = layer_i.evaluateSensitivity(S, xs{i}, options);
    % save sensitivity at layer i for refinement
    if storeSensitivty
        layer_i.sensitivity = S;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
