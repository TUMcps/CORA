function obj = readNNetNetwork(file_path, actFun)
% readNNetNetwork - reads and converts a network saved in nnet format
%
% Syntax:
%    res = neuralNetwork.readNNetNetwork(file_path)
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

% remove comments
lines = lines(~cellfun(@(x) startsWith(x,'//'),lines));

% get number of layers
ind = find(lines{1} == ',');
nrLayers = str2double(lines{1}(1:ind));

% get number of inputs
ind = find(lines{2} == ',');
nrInputs = str2double(lines{2}(1:ind));

% get number of neurons in each layer
nrNeurons = zeros(nrLayers,1);
text = strtrim(lines{2}(ind+1:end));

for i = 1:nrLayers
    ind = find(text == ',');
    nrNeurons(i) = str2double(text(1:ind-1));
    text = strtrim(text(ind+1:end));
end

% loop over all layers and read the weights and biases
W = cell(nrLayers,1); b = cell(nrLayers,1); 
lines = lines(8:end); inpSize = nrInputs;

for i = 1:nrLayers
    % read weight matrix
    W{i} = zeros(nrNeurons(i),inpSize);
    
    for j = 1:nrNeurons(i)
        evalc(['temp = [',lines{j},'];']);
        W{i}(j,:) = temp;
    end
    lines = lines(nrNeurons(i)+1:end);
    
    % read bias
    b{i} = zeros(nrNeurons(i),1);
    
    for j = 1:nrNeurons(i)
        b{i}(j) = str2double(lines{j});
    end
    
    lines = lines(nrNeurons(i)+1:end);
    inpSize = nrNeurons(i);
end

% construct neural network
obj = neuralNetwork.getFromCellArray(W, b, actFun);

% ------------------------------ END OF CODE ------------------------------
