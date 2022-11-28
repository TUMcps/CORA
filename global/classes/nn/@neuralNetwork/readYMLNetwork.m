function obj = readYMLNetwork(file_path)
% readYMLNetwork - reads and converts a network saved in yml format
%
% Syntax:
%    res = neuralNetwork.readYMLNetwork(file_path)
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

% first convert to old neuralNetworkOld
old_nn = neuralnetwork2cora(file_path);

% then convert to new neuralNetwork
obj = neuralNetwork.getFromOldNeuralNetwork(old_nn);

%------------- END OF CODE --------------