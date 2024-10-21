function [gc,gG] = backpropZonotopeBatch_(nn,gc,gG,options,idxLayer)
% backpropZonotopeBatch_ - compute the backpropagation for the previous input
%    with batches of zonotopes without validating the input arguments.
%
% Syntax:
%    [gc, gG] = nn.backpropZonotopeBatch_(gc, gG, options, idxLayer)
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
% See also: neuralNetwork/backprop, neuralNetwork/evaluateZonotopeBatch

% Authors:       Lukas Koller
% Written:       03-August-2023
% Last update:   22-February-2024 (merged options.nn, moved input storage handling from layer to network)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for i = flip(idxLayer)
    layeri = nn.layers{i};
    % Retrieve stored input
    if options.nn.train.backprop
        c = layeri.backprop.store.inc;
        G = layeri.backprop.store.inG;
    end
    [gc,gG] = layeri.backpropZonotopeBatch(c,G,gc,gG,options);
end

end

% ------------------------------ END OF CODE ------------------------------
