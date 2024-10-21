function [nnCopy] = copyNeuralNetwork(obj)
% copyNeuralNetwork - Create a copy of the neural network.
%
% Syntax:
%    [nnCopy] = copyNeuralNetwork(obj)
%
% Inputs:
%    obj - neural network
%
% Outputs:
%    nnCopy - copy of the neural network obj
%    
% References:---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Lukas Koller
% Written:       21-June-2023
% Last update:   ---
% Last revision: ---    

% ------------------------------ BEGIN CODE -------------------------------

% copy layers of the neural network
K = length(obj.layers);
layers = cell(K,1);
for k=1:K
    layers{k} = copy(obj.layers{k});
end
nnCopy = neuralNetwork(layers);
end

% ------------------------------ END OF CODE ------------------------------
