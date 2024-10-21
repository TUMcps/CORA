function [c,G] = evaluateZonotopeBatch(nn,c,G,varargin)
% evaluateZonotopeBatch - evaluate neural network for a batch of zonotopes
%
% Syntax:
%    % prepare network for batch evaluation
%    nn.prepareForZonoBatchEval(c, numInitGens, useApproxError, idxLayer);
%    % compute output sets for entire batch of input sets
%    [c, G] = nn.evaluateZonotopeBatch(c, G);
%
% Inputs:
%    c, G - batch of zonotope; [n,q+1,b] = size([c G]),
%       where n is the number of dims, q the number of generators, and b the batch size
%    options - parameter for neural network evaluation
%    idxLayer - indices of layers that should be evaluated
%
% Outputs:
%    c, G - batch of output sets
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate, neuralNetwork/prepareForZonoBatchEval,
% neuralNetwork/evaluateZonotopeBatch_

% Authors:       Lukas Koller
% Written:       02-August-2023
% Last update:   08-August-2023 (moved code to layers)
%                22-February-2024 (merged options.nn, moved input storage handling from layer to network)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% validate parameters
[options, idxLayer] = setDefaultValues( ...
    {struct, 1:length(nn.layers)}, varargin);
% Set default evaluation parameters.
options = nnHelper.validateNNoptions(options);
% Propagate sets forward through the network.
[c,G] = nn.evaluateZonotopeBatch_(c,G,options,idxLayer);

end

% ------------------------------ END OF CODE ------------------------------
