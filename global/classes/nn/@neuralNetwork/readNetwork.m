function obj = readNetwork(file_path)
% readNetwork - reads and converts a network according to file ending
%
% Syntax:
%    res = neuralNetwork.readNetwork(file_path)
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

[~, ~, type] = fileparts(file_path);

if strcmp(type, '.onnx')
    obj = neuralNetwork.readONNXNetwork(file_path);
elseif strcmp(type, '.nnet')
    obj = neuralNetwork.readNNetNetwork(file_path);
elseif strcmp(type, '.yml')
    obj = neuralNetwork.readYMLNetwork(file_path);
elseif strcmp(type, '.sherlock')
    obj = neuralNetwork.readSherlockNetwork(file_path);
else
    throw(CORAerror('CORA:wrongValue','first',...
        'has to end in ''.onnx'', ''.nnet'', ''.yml'' or ''.sherlock''.'));
end

%------------- END OF CODE --------------