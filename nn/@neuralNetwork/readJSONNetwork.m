function nn = readJSONNetwork(file_path)
% readJSONNetwork - reads and converts a network saved in json format
%
% Syntax:
%    nn = neuralNetwork.readJSONNetwork(file_path)
%
% Inputs:
%    file_path - path to file
% 
% Outputs:
%    nn - neuralNetwork
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/exportAsJSON, neuralNetwork/importFromJSON

% Authors:       Tobias Ladner
% Written:       04-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read json string
jsonstr = fileread(file_path);

% convert to neural network
nn = neuralNetwork.importFromJSON(jsonstr);

% ------------------------------ END OF CODE ------------------------------
