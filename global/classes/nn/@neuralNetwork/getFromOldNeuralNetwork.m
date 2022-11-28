function obj = getFromOldNeuralNetwork(nn_old)
% getFromOldNeuralNetwork - transforms neuralNetworkOld to the new
%    layer-based model.
%
% Syntax:
%    res = neuralNetwork.getFromOldNeuralNetwork(nn_old)
%
% Inputs:
%    nn_old - network of type neuralNetworkOld
%
% Outputs:
%    res - generated network
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork (layer-based), neuralNetworkOld

% Author:       Tobias Ladner
% Written:      28-March-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

layers = cell(length(nn_old.W), 1);
for i = 1:length(nn_old.W)
    % add layers
    layers{2*i-1} = nnLinearLayer(nn_old.W{i}, nn_old.b{i});

    if nn_old.actFun{i} == "sigmoid"
        activation_layer = nnSigmoidLayer();
    elseif nn_old.actFun{i} == "tanh"
        activation_layer = nnTanhLayer();
    elseif nn_old.actFun{i} == "ReLU"
        activation_layer = nnReLULayer();
    elseif nn_old.actFun{i} == "identity"
        activation_layer = nnIdentityLayer();
    else
        throw(CORAerror('CORA:wrongFieldValue','nn_old.actFun',...
            {'sigmoid','tanh','ReLU','identity'}));
    end
    layers{2*i} = activation_layer;
end

obj = neuralNetwork(layers);

%------------- END OF CODE --------------