function [gc,gG] = backpropZonotopeBatch(nn,gc,gG,varargin)
% backpropZonotopeBatch - compute the backpropagation for the previous input
%    with batches of zonotopes
%
% Syntax:
%    % prepare network for batch evaluation
%    nn.prepareForZonoBatchEval(c, numInitGens, useApproxError, idxLayer);
%    % compute output sets for entire batch of input sets
%    [yc, yG] = nn.evaluateZonotopeBatch(c, G);
%    % compute loss and gradients
%    [L, gc, gG] = computeLoss(yc,yG,options);
%    % backpropagate gradients
%    [gc, gG] = nn.backpropZonotopeBatch(gc, gG);
%
% Inputs:
%    gc, gG - batch of zonotope gradients; [n,q+1,b] = size([gc gG]),
%       where n is the number of dims, q the number of generators, and b the batch size
%    options - training parameters
%    idxLayer - indices of layers that should be evaluated
%
% Outputs:
%    gc, gG - zonotope gradients w.r.t the input
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/backprop, neuralNetwork/evaluateZonotopeBatch,
% neuralNetwork/backpropZonotopeBatch_

% Authors:       Lukas Koller
% Written:       03-August-2023
% Last update:   22-February-2024 (merged options.nn, moved input storage handling from layer to network)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% validate parameters
[options, idxLayer] = setDefaultValues( ...
    {struct, 1:length(nn.layers)}, varargin);
% Set default evaluation parameters.
options = nnHelper.validateNNoptions(options,true);
% backpropagte gradients through the network.
[gc,gG] = nn.backpropZonotopeBatch_(gc,gG,options,idxLayer);

end

% ------------------------------ END OF CODE ------------------------------
