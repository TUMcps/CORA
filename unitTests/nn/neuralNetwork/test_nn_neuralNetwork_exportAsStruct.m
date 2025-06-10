function res = test_nn_neuralNetwork_exportAsStruct()
% test_nn_neuralNetwork_exportAsStruct - unit test function for 
%     neuralNetwork/exportAsStruct
%
% Syntax:
%    res = test_nn_neuralNetwork_exportAsStruct()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       10-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% create network
nn = neuralNetwork({ ...
    nnLinearLayer([2 3; 4 5]); ...
    nnSigmoidLayer(); ...
    nnLinearLayer([-1 5; 2 -3]); ...
    nnSigmoidLayer(); ...
});

nnStruct = nn.exportAsStruct();
resvec(end+1) = length(nnStruct) == length(nn.layers);

% empty case
resvec(end+1) = isempty(neuralNetwork().exportAsStruct);

% random case
nn = neuralNetwork.generateRandom();
nnStruct = nn.exportAsStruct();
resvec(end+1) = length(nnStruct) == length(nn.layers);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
