function normWeights(nn,varargin)
% normWeights - normalize the learnable weights of a neural network with 
%   a Lipschitz constant; i.e., call layer.normWeights for all 
%   nnLipConstrLinearLayer
%
% Syntax:
%    nn.normWeights(x,idxLayer)
%
% Inputs:
%    nn - object of class neuralNetwork
%    options - 
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
% See also: neuralNetwork, nnLipConstrLinearLayer

% Authors:       Lukas Koller
% Written:       04-December-2023
% Last update:   18-August-2025 (enumerate layers)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Enumerate all layers.
[layers,~] = nn.enumerateLayers();

% Validate parameters.
[options,idxLayer] = setDefaultValues({struct,1:length(layers)}, varargin);
options = nnHelper.validateNNoptions(options);
   
for i=idxLayer
    % Obtain the i-th layer.
    layeri = layers{i};
    if isa(layeri,'nnLipConstrLinearLayer') ...
            && ~isempty(layeri.getLearnableParamNames())
        % Normalize the weights of the layer.
        layeri.normWeights(options);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
