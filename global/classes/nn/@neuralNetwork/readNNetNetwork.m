function obj = readNNetNetwork(file_path)
% readNNetNetwork - reads and converts a network saved in nnet format
%
% Syntax:
%    res = neuralNetwork.readNNetNetwork(file_path)
%
% Inputs:
%    file_path: path to file
% 
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralnetwork2cora

% Author:       Tobias Ladner
% Written:      30-March-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% TODO: remove dependency on neuralNetworkOld

% first convert to old neuralNetwork
old_nn = neuralnetwork2cora(file_path);

% then convert to new neuralNetwork
obj = neuralNetwork.getFromOldNeuralNetwork(old_nn);

%------------- END OF CODE --------------