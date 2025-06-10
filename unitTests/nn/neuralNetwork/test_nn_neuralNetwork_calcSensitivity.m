function [res] = test_nn_neuralNetwork_calcSensitivity()
% test_nn_neuralNetwork_calcSensitivity - tests the calcSensitivity function 
%    
%
% Syntax:
%    res = test_nn_neuralNetwork_calcSensitivity()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       23-April-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Specify number of input and output dimensions.
n0 = 5;
nK = 7;
% Generate a random neural network.
nn = neuralNetwork.generateRandom( ...
    'NrInputs',n0, ...
    'NrOutputs',nK, ...
    'NrLayers',3, ...
    'NrHiddenNeurons',17 ...
);
% Append a softmax layer.
nn.layers{end+1} = nnSoftmaxLayer;

% Specify a batch size.
bSz = 13;
% Generate a random input.
x = rand([n0 bSz]);
% Compute the output.
y = nn.evaluate(x);

% Calculate the sensitivity.
S = nn.calcSensitivity(x,struct,true);

% Check the dimensions of the sensitivity matrix.
assert(all(size(S) == [nK n0 bSz]));

% Generate a second random input.
x_ = rand([n0 bSz]);
% Compute the difference between the two inputs.
dx = x - x_;
% Compute the new output.
y_ = nn.evaluate(x + dx);
% Compute expected difference based on the sensitivity.
dy = pagemtimes(S,permute(dx,[1 3 2]));
% Check if the directions of the sensitivity matrix are correct.
assert(all(sign(y + dy) == sign(y_),'all'));

% Extract the sensitivity matrix of the last layer.
Sk = nn.layers{end-1}.sensitivity;
% Check the dimensions of the sensitivity matrix.
assert(all(size(Sk) == [nK nK bSz]));

res = 1;

end

% ------------------------------ END OF CODE ------------------------------
