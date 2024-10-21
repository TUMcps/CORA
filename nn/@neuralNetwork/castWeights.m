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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% validate parameters
[idxLayer] = setDefaultValues({1:length(nn.layers)}, varargin);
   
for i = idxLayer
    layeri = nn.layers{i};
    % move all learnable parameters to gpu
    names = layeri.getLearnableParamNames();
    for j=1:length(names)
        % cast learnable weights
        layeri.(names{j}) = cast(layeri.(names{j}),'like',x);
    end
    % Notify layer
    layeri.castWeights(x);
end

end

% ------------------------------ END OF CODE ------------------------------
