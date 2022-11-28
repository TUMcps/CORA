function obj = generateRandom(varargin)
% generateRandom - creates a random layer-based network
%
% Syntax:
%    obj = neuralNetwork.generateRandom()
%    obj = neuralNetwork.generateRandom('NrInputs',nrOfInputs)
%    obj = neuralNetwork.generateRandom('NrInputs',nrOfInputs,...
%       'NrOutputs',nrOfOutputs)
%    obj = neuralNetwork.generateRandom('NrInputs',nrOfInputs,...
%       'NrOutputs',nrOfOutputs,'ActivationFun',actFun)
%    obj = neuralNetwork.generateRandom('NrInputs',nrOfInputs,...
%       'NrOutputs',nrOfOutputs,'ActivationFun',actFun,'NrLayers',nrLayers)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'NrInputs',nrOfInputs> - number of inputs
%       <'NrOutputs',nrOfOutputs> - number of outputs
%       <'ActivationFun',actFun> - type of activation functions
%           actFun has to be either {'ReLU', 'sigmoid', 'tanh'}
%       <'NrLayers',nrLayers> - number of layers
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Author:       Tobias Ladner
% Written:      28-March-2022
% Last update:  28-November-2022 (name-value pair syntax)
% Last revision:---

%------------- BEGIN CODE --------------

% validate parameters
% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'NrInputs','NrOutputs','NrLayers','ActivationFun'});
    % number of inputs given?
    [NVpairs,nrInputs] = readNameValuePair(NVpairs,'NrInputs');
    % number of outputs given?
    [NVpairs,nrOutputs] = readNameValuePair(NVpairs,'NrOutputs');
    % activation function given?
    [NVpairs,actFun] = readNameValuePair(NVpairs,'ActivationFun');
    % number of layers given?
    [NVpairs,nrLayers] = readNameValuePair(NVpairs,'NrLayers');
end

% old
if isempty(nrInputs)
    nrInputs = randi([1, 5]);
end
if isempty(nrOutputs)
    nrOutputs = randi([1, 5]);
end
if isempty(actFun)
    actFun = "sigmoid";
end
if isempty(nrLayers)
    nrLayers = randi([1, 5]);
end

% determine neurons in each layer
neurons = zeros(1, 1+nrLayers);
neurons(1) = nrInputs;
for i = 1:nrLayers - 1
    neurons(1+i) = randi([1, 20]);
end
neurons(end) = nrOutputs;

% create layers
layers = cell(2*(length(neurons) - 1), 1);

scale = 1;
for i = 1:length(neurons) - 1
    % add linear layer
    W = rand(neurons(i+1), neurons(i)) * scale - scale / 2;
    b = rand(neurons(i+1), 1) * scale - scale / 2;
    layers{2*i-1} = nnLinearLayer(W, b);

    % add activation layers
    layer = nnActivationLayer.instantiateFromString(actFun);
    layers{2*i} = layer;
end

% create neuralNetwork
obj = neuralNetwork(layers);

%------------- END OF CODE --------------