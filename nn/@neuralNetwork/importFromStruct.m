function nn = importFromStruct(nnStruct)
% importFromStruct - imports a network from struct
%
% Syntax:
%    nn = neuralNetwork.importFromStruct(nnStruct)
%
% Inputs:
%    nnStruct - struct
%
% Outputs:
%    nn - neuralNetwork
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/exportAsStruct

% Authors:       Tobias Ladner
% Written:       10-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read layers
layers = cell(length(nnStruct),1);
for i = 1:length(nnStruct)
    layers{i} = nnLayer.importFromStruct(nnStruct(i));
end

% init neuralNetwork
nn =  neuralNetwork(layers);

end

% ------------------------------ END OF CODE ------------------------------
