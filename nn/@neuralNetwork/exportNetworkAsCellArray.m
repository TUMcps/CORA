function res = exportNetworkAsCellArray(nn, file_path, varargin)
% exportNetworkAsCellArray - exports the given neural network as
%    cell array containing weights, bias, and activation functions
%
% Syntax:
%    res = exportNetworkAsCellArray(nn,file_path)
%
% Inputs:
%    nn - neuralNetwork
%    file_path - path to file
%
% Outputs:
%    res - logical
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/getFromCellArray, neuralNetwork/exportONNXNetwork

% Authors:       Tobias Ladner
% Written:       02-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 2
    file_path = './network.mat';
end

% normalize network
nn = nn.getNormalForm();
K = numel(nn.layers);

% collect properties
W = cellfun(@(layer) layer.W, nn.layers(1:2:K), 'UniformOutput', false);
b = cellfun(@(layer) layer.b, nn.layers(1:2:K), 'UniformOutput', false);
actFun = cellfun(@(layer) lower(strrep(strrep(class(layer),'nn',''),'Layer','')), nn.layers(2:2:K), 'UniformOutput', false);

% save properties
save(file_path, 'W', 'b', 'actFun')

% result
res = true;

% ------------------------------ END OF CODE ------------------------------
