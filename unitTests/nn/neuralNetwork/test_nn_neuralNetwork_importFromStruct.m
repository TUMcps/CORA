function res = test_nn_neuralNetwork_importFromStruct()
% test_nn_neuralNetwork_importFromStruct - unit test function for 
%     neuralNetwork/importAsStruct
%
% Syntax:
%    res = test_nn_neuralNetwork_importFromStruct()
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
nn2 = neuralNetwork.importFromStruct(nnStruct);

x = ones(2,1);
resvec(end+1) = all(withinTol(nn.evaluate(x), nn2.evaluate(x),1e-6));

% empty case
nn = neuralNetwork();
nnStruct = nn.exportAsStruct();
nn2 = neuralNetwork.importFromStruct(nnStruct);

x = ones(2,1);
resvec(end+1) = all(withinTol(nn.evaluate(x), nn2.evaluate(x),1e-6));

% random case
nn = neuralNetwork.generateRandom();
nnStruct = nn.exportAsStruct();
nn2 = neuralNetwork.importFromStruct(nnStruct);

x = ones(nn.neurons_in,1);
resvec(end+1) = all(withinTol(nn.evaluate(x), nn2.evaluate(x),1e-6));

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
