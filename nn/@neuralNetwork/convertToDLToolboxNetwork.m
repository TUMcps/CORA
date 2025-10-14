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
%    nn_dlt - dlnetwork
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       30-April-2024
% Last update:   09-May-2025 (TL, added conv + batch norm + avg pool)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init layers
if isempty(nn.neurons_in)
    throw(CORAerror('CORA:converterIssue','Unable to determine number of input neurons. See also neuralNetwork/setInputSize.'))
end
layers_cora = nn.layers;

% init input of dlt network
if any(cellfun(@(layer_cora) isa(layer_cora,'nnConv2DLayer'), layers_cora))
    % image input
    layers_dlt = inputLayer(layers_cora{1}.inputSize,'SSC');
else
    % feature input
    layers_dlt = featureInputLayer(nn.neurons_in);
end

% iterate over all layers
for k=1:numel(layers_cora)
    layer_cora_k = layers_cora{k};

    switch class(layer_cora_k)

        % linear ---

        case 'nnLinearLayer'
            W = layer_cora_k.W;
            b = layer_cora_k.b;

            % check approximation error
            if ~representsa(layer_cora_k.d,'emptySet')
                CORAwarning('CORA:nn','Converting linear layer with non-empty approximation error. Adding center to bias...')
                b = b + center(layer_cora_k.d);
            end

            % init dlt layer
            layer_dlt_k = fullyConnectedLayer(size(W,1));
            layer_dlt_k.Weights = W;
            layer_dlt_k.Bias = b;

        case 'nnElementwiseAffineLayer'
            scale = layer_cora_k.scale;
            offset = layer_cora_k.offset;
            layer_dlt_k = nnet.onnx.layer.ElementwiseAffineLayer('', scale, offset);

        % convolutional ---

        case 'nnConv2DLayer'
            % conv 2d
            filterSize = size(layer_cora_k.W,1:2);
            numFilters = size(layer_cora_k.W,4);
            layer_dlt_k = convolution2dLayer(filterSize,numFilters,...
                'Weights',layer_cora_k.W, ...
                'Bias',reshape(layer_cora_k.b,1,1,[]), ...
                'Stride',layer_cora_k.stride, ...
                'Padding',layer_cora_k.padding, ...
                'DilationFactor',layer_cora_k.dilation ...
                );
            
        case 'nnBatchNormLayer'
            % batch norm
            layer_dlt_k = batchNormalizationLayer;
            layer_dlt_k.Scale = layer_cora_k.scale;
            layer_dlt_k.Offset = layer_cora_k.offset;
            layer_dlt_k.TrainedMean = layer_cora_k.movMean;
            layer_dlt_k.TrainedVariance = layer_cora_k.movVar;
            layer_dlt_k.Epsilon = layer_cora_k.epsilon;

        case 'nnAvgPool2DLayer'
            % avg pooling
            layer_dlt_k = averagePooling2dLayer( ...
                layer_cora_k.poolSize, ...
                'Stride',layer_cora_k.stride ...
            );

        case 'nnReshapeLayer'
            % assuming flattening; requires adaptation for other cases
            layer_dlt_k = flattenLayer;

        % nonlinear/activation ---

        case 'nnReLULayer'
            layer_dlt_k = reluLayer;

        case 'nnTanhLayer'
            layer_dlt_k = tanhLayer;

        case 'nnSigmoidLayer'
            layer_dlt_k = sigmoidLayer;

        otherwise % throw error ---
            throw(CORAerror('CORA:converterIssue', sprintf("Layer conversion of '%s' not implemented yet.", class(layer_cora_k))))

    end

    % add name
    layer_dlt_k.Name = layer_cora_k.name;

    % add to dlt layers list
    layers_dlt(end+1) = layer_dlt_k;
end

% init DLT network
nn_dlt = dlnetwork(layers_dlt);

end


% ------------------------------ END OF CODE ------------------------------
