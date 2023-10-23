function obj = readNetwork(file_path,varargin)
% readNetwork - reads and converts a network according to file ending
%
% Syntax:
%    res = neuralNetwork.readNetwork(file_path)
%
% Inputs:
%    file_path: path to file
%    varargin: further input parameter
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       30-March-2022
% Last update:   30-November-2022 (inputArgsCheck)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% validate input
if ~isfile(file_path)
    % test if file is already on path
    found_file = which(file_path);
    if isempty(found_file)
        throw(CORAerror('CORA:fileNotFound', file_path));
    else
        file_path = found_file;
    end
end
[~, ~, ext] = fileparts(file_path);
possibleExtensions = {'.onnx', '.nnet', '.yml', '.sherlock'};
inputArgsCheck({{ext, 'str', possibleExtensions}});

% redirect to specific read function
if strcmp(ext, '.onnx')
    obj = neuralNetwork.readONNXNetwork(file_path,varargin{:});
elseif strcmp(ext, '.nnet')
    obj = neuralNetwork.readNNetNetwork(file_path,varargin{:});
elseif strcmp(ext, '.yml')
    obj = neuralNetwork.readYMLNetwork(file_path,varargin{:});
elseif strcmp(ext, '.sherlock')
    obj = neuralNetwork.readSherlockNetwork(file_path,varargin{:});
else
    throw(CORAerror('CORA:wrongValue','first',...
        sprintf('has to end in %s', strjoin(possibleExtensions, ', ')) ...
    ));
end

% ------------------------------ END OF CODE ------------------------------
