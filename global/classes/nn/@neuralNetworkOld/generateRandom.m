function obj = generateRandom(varargin)
% generateRandom - generates a random neural network
%
% Syntax:  
%    obj = generateRandom()
%    obj = generateRandom('NrInputs',nrOfInputs)
%    obj = generateRandom('NrInputs',nrOfInputs,'NrOutputs',nrOfOutputs)
%    obj = generateRandom('NrInputs',nrOfInputs,'NrOutputs',nrOfOutputs,...
%       'ActivationFun',actFun)
%    obj = generateRandom('NrInputs',nrOfInputs,'NrOutputs',nrOfOutputs,...
%       'ActivationFun',actFun,'NrLayers',layers)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'NrInputs',nrOfInputs> - number of inputs
%       <'NrOutputs',nrOfOutputs> - number of outputs
%       <'ActivationFun',actFun> - type of activation functions
%           actFun has to be {'ReLU', 'sigmoid', 'tanh'}
%       <'NrLayers',layers> - number of layers
%
% Outputs:
%    obj - generated neuralNetworkOld object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetworkOld

% Author:       Niklas Kochdumper
% Written:      17-September-2021             
% Last update:  19-May-2022 (MW, name-value pair syntax)
% Last revision:---

%------------- BEGIN CODE --------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'NrInputs','NrOutputs','ActivationFun','NrLayers'});
    % number of inputs given?
    [NVpairs,nrInputs] = readNameValuePair(NVpairs,'NrInputs');
    % number of outputs given?
    [NVpairs,nrOutputs] = readNameValuePair(NVpairs,'NrOutputs');
    % activation function given?
    [NVpairs,actFun] = readNameValuePair(NVpairs,'ActivationFun');
    % number of layers given?
    [NVpairs,layers] = readNameValuePair(NVpairs,'NrLayers');
end

% default value for number of inputs
if isempty(nrInputs)
    nrInputs = randi([1,20]);
end

% default value for number of outputs
if isempty(nrOutputs)
    nrOutputs = randi([1,20]);
end

% default activation function
validActFuns = {'ReLU','sigmoid','tanh'};
if isempty(actFun)
    actFun = validActFuns{randi([1,3])};
elseif ~ismember(actFun,validActFuns)
    throw(CORAerror('CORA:wrongValue','ActivationFun',...
        {'ReLU','sigmoid','tanh'}));
end

% default number of layers
if isempty(layers)
    layers = randi([1,10]);
end

% set number of neurons for all layers
nrNeurons = zeros(layers+1,1);
nrNeurons(1) = nrInputs;
nrNeurons(end) = nrOutputs;

for i = 2:layers
    nrNeurons(i) = randi([1,100]);
end

% set values for weights and constant offsets
W = cell(layers,1); b = cell(layers,1);

for i = 1:layers
    W{i} = 0.5 - rand(nrNeurons(i+1),nrNeurons(i));
    b{i} = 0.5 - rand(nrNeurons(i+1),1);
end

% construct neuralNetworkOld object
obj = neuralNetworkOld(W,b,actFun);

%------------- END OF CODE --------------