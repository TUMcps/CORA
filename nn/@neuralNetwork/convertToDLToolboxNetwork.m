function nn_dlt = convertToDLToolboxNetwork(nn)
% convertToDLToolboxNetwork - converts a CORA network to a network
%    from the the Deep Learning Toolbox
%
% Syntax:
%    res = convertToDLToolboxNetwork(nn)
%
% Inputs:
%    nn - dlnetwork
%
% Outputs:
%    nn_dlt - 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       30-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init layers
if isempty(nn.neurons_in)
    throw(CORAerror('CORA:converterIssue','Unable to determine number of input neurons. See also neuralNetwork/setInputSize.'))
end
layers_dlt = featureInputLayer(nn.neurons_in);
layers_cora = nn.layers;

% iterate over all layers
for k=1:numel(layers_cora)
    layer_cora_k = layers_cora{k};

    switch class(layer_cora_k)

        case 'nnLinearLayer'
            % read properties
            W = layer_cora_k.W;
            b = layer_cora_k.b;

            % init dlt layer
            layer_dlt_k = fullyConnectedLayer(size(W,1));
            layer_dlt_k.Weights = W;
            layer_dlt_k.Bias = b;

        case 'nnReLULayer'
            layer_dlt_k = reluLayer;

        case 'nnTanhLayer'
            layer_dlt_k = tanhLayer;

        case 'nnSigmoidLayer'
            layer_dlt_k = sigmoidLayer;

        case 'nnElementwiseAffineLayer'
            scale = layer_cora_k.scale;
            offset = layer_cora_k.offset;
            layer_dlt_k = nnet.onnx.layer.ElementwiseAffineLayer('', scale, offset);

        otherwise
            throw(CORAerror('CORA:converterIssue', sprintf("Layer conversion of '%s' not implemented yet.", class(layer_cora_k))))

    end

    % add name
    layer_dlt_k.Name = layer_cora_k.name;

    % add to dlt layers list
    layers_dlt(end+1) = layer_dlt_k;
end


% init DLT network
nn_dlt = dlnetwork(layers_dlt);


% ------------------------------ END OF CODE ------------------------------
