function obj = readSherlockNetwork(file_path, actFun)
% readSherlockNetwork - reads and converts a network saved in sherlock
%    format
%
% Syntax:
%    res = neuralNetwork.readSherlockNetwork(file_path)
%
% Inputs:
%    file_path - path to file
%    actFun - string of activation function
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/getFromCellArray

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       12-November-2021
% Last update:   30-March-2022
%                30-November-2022 (removed neuralNetworkOld)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% actFun will be checked in neuralNetwork/getFromCellArray at the end
if nargin < 2
    % support legacy
    actFun = 'ReLU';
end

% read text from file
text = fileread(file_path);
lines = strsplit(text,'\n');

% get network properties
nrInputs = str2double(strtrim(lines{1}));
nrOutputs = str2double(strtrim(lines{2}));
hiddenLayers = str2double(strtrim(lines{3}));

% get number of neurons in each layer
nrNeurons = cell(hiddenLayers + 2,1);
nrNeurons{1} = nrInputs;
nrNeurons{end} = nrOutputs;

for i = 1:hiddenLayers
    nrNeurons{i+1} = str2double(strtrim(lines{3+i}));
end

% initialization
cnt = 3 + hiddenLayers;
W = cell(hiddenLayers+1,1);
b = cell(hiddenLayers+1,1);

% loop over all layers
for i = 1:length(nrNeurons)-1
    
    % initialization
    temp = zeros(nrNeurons{i+1},nrNeurons{i}+1);
   
    % read data
    for k = 1:nrNeurons{i+1}
        offset = (k-1)*(nrNeurons{i}+1);
        for j = 1:nrNeurons{i}+1
            temp(k,j) = str2double(strtrim(lines{cnt+offset+j}));
        end
    end
    cnt = cnt + (nrNeurons{i}+1)*nrNeurons{i+1};
    
    % get weight matrix and bias vector
    W{i} = temp(:,1:end-1);
    b{i} = temp(:,end);
end

% construct neural network
obj = neuralNetwork.getFromCellArray(W,b,actFun);

% ------------------------------ END OF CODE ------------------------------
