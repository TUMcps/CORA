function nn = importFromJSON(jsonstr)
% importFromJSON - imports a network from json
%
% Syntax:
%    nn = neuralNetwork.importFromJSON(jsonstr)
%
% Inputs:
%    jsonstr - json string
%
% Outputs:
%    nn - neuralNetwork
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/exportAsJSON

% Authors:       Tobias Ladner
% Written:       10-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% decode json
nnStruct = jsondecode(jsonstr);

% import network
nn = neuralNetwork.importFromStruct(nnStruct);

end

% ------------------------------ END OF CODE ------------------------------
