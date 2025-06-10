function res = test_nn_neuralNetwork_readGNNnetwork()
% test_nn_neuralNetwork_readGNNnetwork - unit test function for reading a
%     GNN netowrk and data
%
% Syntax:
%    res = test_nn_neuralNetwork_readGNNnetwork()
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
% See also: neuralNetwork/evaluate

% Authors:       Tobias Ladner
% Written:       27-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

models = {'Enzymes', 'Proteins'};

for m = 1:numel(models)
    model = models{m};
    modelpath = [CORAROOT '/models/Cora/nn/gnn/' model];
    fprintf('Reading network and data in %s..\n', model)
    nn = neuralNetwork.readGNNnetwork([modelpath filesep 'model_export.json']);
    data = neuralNetwork.readGNNdata([modelpath filesep 'data_export.json']);

    fprintf('Success. Model has %i message passing steps. Also loaded dataset with %i graphs.\n\n', nn.getNumMessagePassingSteps(), height(data));
end

res = true;

end

% ------------------------------ END OF CODE ------------------------------
