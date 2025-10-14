function castWeights(nn,x,varargin)
% castWeights - cast learnable weights of a neural network
%
% Syntax:
%    nn.castWeights(x,idxLayer)
%
% Inputs:
%    nn - object of class neuralNetwork
%    x - instance of data type
%    idxLayer - indices of layers that should be evaluated
%
% Outputs:
%    -
% 
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Lukas Koller
% Written:       04-December-2023
% Last update:   18-August-2025 (enumerate layers)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Enumerate all layers.
[layers,~] = nn.enumerateLayers();

% Validate parameters.
[idxLayer] = setDefaultValues({1:length(layers)}, varargin);
   
for i=idxLayer
    layeri = layers{i};
    % move all learnable parameters to gpu
    names = layeri.getParamNames();
    for j=1:length(names)
        % cast learnable weights
        layeri.(names{j}) = cast(layeri.(names{j}),'like',x);
    end
    % Notify layer
    layeri.castWeights(x);
end

end

% ------------------------------ END OF CODE ------------------------------
