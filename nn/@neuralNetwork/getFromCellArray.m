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
%   W = {[ 0.837 -0.893 ; 0.475 -0.895 ; 1.252 0.313 ],[ 0.667 -0.236 0.885 ; 0.828 0.653 0.085 ; 0.007 1.965 -0.575 ], [ 0.500 0.778 ; -0.484 0.924 ; 0.238 0.588 ]};
%   b = {[ 1.378 ; 1.851 ; -1.898 ], [ 1.378 ; 1.851 ; -1.898 ], [ -0.357 ; -0.336 ; 0.250 ]};
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
