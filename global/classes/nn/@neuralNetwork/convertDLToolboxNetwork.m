function obj = convertDLToolboxNetwork(dltoolbox_layers, verbose)
% convertDLToolboxNetwork - converts a network from the Deep Learning
%    Toolbox to a CORA neuralNetwork for verification
%
% Syntax:
%    res = neuralNetwork.convertDLToolboxNetwork(dltoolbox_layers)
%    res = neuralNetwork.convertDLToolboxNetwork(dltoolbox_layers,verbose)
%
% Inputs:
%    dltoolbox_layers - layer array (e.g. dltoolbox_nn.Layers)
%    verbose - true/false whether information should be displayed
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralnetwork2cora

% Author:       Tobias Ladner
% Written:      30-March-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

n = size(dltoolbox_layers, 1);

m = 0;
layers = cell(n, 1);

if verbose
    disp("Converting Deep Learning Toolbox Model to neuralNetwork...")
end

for i = 1:n
    dlt_layer = dltoolbox_layers(i);
    if verbose
        fprintf("#%d: %s\n", i, class(dlt_layer))
    end

    % handle different types of layers
    if isa(dlt_layer, 'nnet.cnn.layer.FullyConnectedLayer')
        m = m + 1;

        W = double(dlt_layer.Weights);
        b = double(dlt_layer.Bias);
        layers{m} = nnLinearLayer(W, b, dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.onnx.layer.ElementwiseAffineLayer')
        s = double(dlt_layer.Scale);
        o = squeeze((dlt_layer.Offset));

        if size(o, 1) >= 1 && size(o, 2) == 1
            m = m + 1;
            layers{m} = nnElementwiseAffineLayer(s, o, dlt_layer.Name);
        elseif size(o, 1) == 1 && size(o, 2) > 1
            m = m + 1;
            layers{m} = nnElementwiseAffineLayer(s, o', dlt_layer.Name);
        else
            throw(CORAerror('CORA:converterIssue',...
                ['Conversion failed due to incompatibility with ', ...
                '"Elementwise Affine Layer"!']));
        end

        % activation functions
    elseif isa(dlt_layer, 'nnet.cnn.layer.ReLULayer')
        m = m + 1;
        layers{m} = nnReLULayer(dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.LeakyReLULayer')
        alpha = double(dlt_layer.Scale);
        m = m + 1;
        layers{m} = nnLeakyReLULayer(alpha, dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.TanhLayer')
        m = m + 1;
        layers{m} = nnTanhLayer(dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.SigmoidLayer')
        m = m + 1;
        layers{m} = nnSigmoidLayer(dlt_layer.Name);

    else
        % show warning
        if verbose
            warning("Skipping '%s'. Not implemented in cora yet!.", class(dlt_layer))
        end
    end
end


layers = layers(1:m, :);
if verbose
    disp("Following layers were read successfully:")
    disp(layers)
    fprintf("(%d layers)\n", size(layers, 1))
end

% instantiate neural network
obj = neuralNetwork(layers);

%------------- END OF CODE --------------