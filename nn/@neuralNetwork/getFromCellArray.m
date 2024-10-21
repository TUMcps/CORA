function obj = getFromCellArray(W, b, actFun)
% getFromCellArray - instantiates a neuralNetwork from the given weight,
%    bias, and activation function stored in cell arrays
%    (legacy function, move to the neuralNetwork/read... function or 
%     instantiate layers directly)
%
% Syntax:
%    res = neuralNetwork.getFromCellArray(W, b, {actFun1, actFun2, ...})
%    res = neuralNetwork.getFromCellArray(W, b, actFun)
%
% Inputs:
%    W - cell array holding the weight matrix per linear layer
%    b - cell array holding the bias vectors per linear layer
%    actFun:
%       - cell array holding the actFun string per nonlinear layer
%       - actFun string for all nonlinear layer
%       (see possible strings in nnActivationLayer/instantiateFromString)
%
% Outputs:
%    res - generated network
%
% Example:
%   W = {rand(3, 2), rand(3, 3), rand(3, 2)};
%   b = {rand(3, 1), rand(3, 1), rand(3, 1)};
%   actFun = {'relu', 'relu', 'identity'};
%   nn = neuralNetwork.getFromCellArray(W, b, actFun);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork, nnActivationLayer/instantiateFromString

% Authors:       Tobias Ladner
% Written:       30-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% pre-process second input possibility
K = length(W);
if ~isa(actFun, 'cell')
    actFunStr = actFun;
    actFun = cell(1, K);
    actFun(:) = {actFunStr};
end

% validate input
inputArgsCheck({ ...
    {W, 'att', {'cell'}}; ...
    {b, 'att', {'cell'}}; ...
    {actFun, 'att', {'cell'}}; ...
})
if ~(length(b) == K)
    throw(CORAerror('CORA:dimensionMismatch', W, b))
end
K = length(W);
if ~(length(b) == K)
    throw(CORAerror('CORA:dimensionMismatch', W, actFun))
end

% instantiate layers
layers = cell(2*K, 1);
for k=1:K
    layers{2*k-1} = nnLinearLayer(W{k}, b{k});
    layers{2*k} = nnActivationLayer.instantiateFromString(actFun{k});
end
obj = neuralNetwork(layers);

% ------------------------------ END OF CODE ------------------------------
