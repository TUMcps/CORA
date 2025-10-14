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
%    idxLayer - indices of layers to be evaluated
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
narginchk(2,5);
[options,storeSensitivty,idxLayer] = setDefaultValues(...
    {struct,true,1:length(obj.layers)}, varargin);
options = nnHelper.validateNNoptions(options);

% Enable backpropagation to store the inputs to each layer.
options.nn.train.backprop = true;

% Do a forward propagation to store the inputs.
y = obj.evaluate(x,options,idxLayer);
% Obtain number of output dimensions and batch size.
[nK,bSz] = size(y);

% Initialize the sensitivity in for the output, i.e., identity matrix.
S = repmat(eye(nK,'like',y),1,1,bSz);

% Enable storing the sensitivity.
options.nn.store_sensitivity = storeSensitivty;

% Back-propagate the sensitivity matrix.
for i=flip(idxLayer)
    % Obtain the i-th layer.
    layeri = obj.layers{i};
    % Evaluate the sensitivity.
    S = layeri.evaluateSensitivity(S,options);
end

end

% ------------------------------ END OF CODE ------------------------------
