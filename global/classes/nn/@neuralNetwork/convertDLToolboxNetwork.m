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
% See also: -

% Authors:       Tobias Ladner, Lukas Koller
% Written:       30-March-2022
% Last update:   05-June-2022 (LK, Conv, Pool)
%                17-January-2023 (TL, Reshape)
%                23-November-2023 (TL, bug fix with scalar element-wise operation)
% Last revision: 17-August-2022
%                25-July-2023 (TL, nnElementwiseAffineLayer)

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 2
    verbose = false;
end

n = size(dltoolbox_layers, 1);

layers = {};
inputSize = [];
currentSize = [];

if verbose
    disp("Converting Deep Learning Toolbox Model to neuralNetwork...")
end

for i = 1:n
    dlt_layer = dltoolbox_layers(i);
    if verbose
        fprintf("#%d: %s\n", i, class(dlt_layer))
    end

    % handle different types of layers
    if isa(dlt_layer, 'nnet.cnn.layer.ImageInputLayer') || ...
            isa(dlt_layer, 'nnet.cnn.layer.FeatureInputLayer')
        inputSize = dlt_layer.InputSize;
        if length(inputSize) == 1
            inputSize = [inputSize, 1];
        elseif length(inputSize) == 3
            % channel dimension should be last: [h,w,c]
            inputSize = sort(inputSize, 'descend');
        end
        currentSize = inputSize;
        continue;

        % Normal Layers ---------------------------------------------------

    elseif isa(dlt_layer, 'nnet.cnn.layer.FullyConnectedLayer')
            W = double(dlt_layer.Weights);
            b = double(dlt_layer.Bias);
            layers{end+1} = nnLinearLayer(W, b, dlt_layer.Name);
    
    elseif isa(dlt_layer, 'nnet.onnx.layer.ElementwiseAffineLayer') 
        s = double(dlt_layer.Scale);
        o = double(dlt_layer.Offset);

        % fix dimensions for [h,w,c] inputs
        if length(currentSize) == 3
            if ~isscalar(s)
                % fix if all values are equal
                if ~isempty(s) && all(s(1) == s,'all')
                    s = s(1);
                % try to fix scaling factor
                elseif length(currentSize) ~= length(size(s))
                    % assuming scaling in channel dimension
                    s = reshape(s, 1, 1, []);
                    s = repmat(s, [currentSize(1:2), 1]);
                elseif ~all(currentSize == size(s))
                    % reshaping to current size
                    s = reshape(s,currentSize);
                end
            end
            if ~isscalar(o)
                % fix if all values are equal
                if ~isempty(o) && all(o(1) == o,'all')
                    o = o(1);
                % try to fix offset vector
                elseif length(currentSize) ~= length(size(o))
                    % assuming scaling in channel dimension
                    o = reshape(o, 1, 1, []);
                    o = repmat(o, [currentSize(1:2), 1]);
                elseif ~all(currentSize == size(o))
                    % reshaping to current size
                    o = reshape(o,currentSize);
                end
            end
        end

        % should be column vector
        s = reshape(s, [], 1);
        o = reshape(o, [], 1);

        layers{end+1} = nnElementwiseAffineLayer(s, o, dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.BatchNormalizationLayer')
        % can be converted to elementwise layers
        % https://onnx.ai/onnx/operators/onnx__BatchNormalization.html

        mean = dlt_layer.TrainedMean;
        var = dlt_layer.TrainedVariance;
        epsilon = dlt_layer.Epsilon;
        scale = dlt_layer.Scale;
        bias = dlt_layer.Offset;

        % (x-mean) / sqrt(var+epsilon) * scale + B
        % = x / sqrt(var+epsilon) * scale + (B - mean / sqrt(var+epsilon) * scale)
        % thus,

        final_scale =  1 ./ sqrt(var+epsilon) .* scale;
        final_offset = (bias - mean .* final_scale);

        layers{end+1} = nnElementwiseAffineLayer(final_scale, final_offset, dlt_layer.Name);
        
    elseif isa(dlt_layer, 'nnet.cnn.layer.AveragePooling2DLayer')
        poolSize = dlt_layer.PoolSize;
        padding = dlt_layer.PaddingSize;
        stride = dlt_layer.Stride;
        dilation = [1, 1];

        layers{end+1} = nnAvgPool2DLayer(poolSize, padding, stride, dilation, dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.Convolution2DLayer')
        W = double(dlt_layer.Weights);
        b = double(reshape(dlt_layer.Bias, [], 1));
        padding = dlt_layer.PaddingSize;
        stride = dlt_layer.Stride;
        dilation = dlt_layer.DilationFactor;
        layers{end+1} = nnConv2DLayer(W, b, padding, stride, dilation, dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.MaxPooling2DLayer')
        poolSize = dlt_layer.PoolSize;
        stride = dlt_layer.Stride;
        layers{end+1} = nnMaxPool2DLayer(poolSize, stride, dlt_layer.Name);

        % activation functions
    elseif isa(dlt_layer, 'nnet.cnn.layer.ReLULayer')
        layers{end+1} = nnReLULayer(dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.LeakyReLULayer')
        alpha = double(dlt_layer.Scale);
        layers{end+1} = nnLeakyReLULayer(alpha, dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.TanhLayer')
        layers{end+1} = nnTanhLayer(dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.SigmoidLayer')
        layers{end+1} = nnSigmoidLayer(dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.onnx.layer.IdentityLayer')
        % ignore
        continue;

    elseif isa(dlt_layer, 'nnet.cnn.layer.RegressionOutputLayer')
        % ignore
        continue;

        % Custom Layers -----------------------------------------------

    elseif contains(lower(class(dlt_layer)), 'flatten') || ...
            contains(lower(class(dlt_layer)), 'reshape')
        % flatten
        
        idx = dlarray(1:prod(currentSize));
        idx = reshape(idx, currentSize);
        idx_out = dlt_layer.predict(idx);
        idx_out = double(extractdata(idx_out));
        layers{end+1} = nnReshapeLayer(idx_out, dlt_layer.Name);

    elseif isa(dlt_layer, 'test_sat.Add_To_Sub1018')
        % for VNN Comp
        W1 = dlt_layer.Operation_1_MatMul_W;
        layers{end+1} = nnLinearLayer(W1, zeros(size(W1, 1)));

        W2 = dlt_layer.Operation_2_MatMul_W;
        layers{end+1} = nnLinearLayer(W2, zeros(size(W2, 1)));

        W3 = dlt_layer.Operation_3_MatMul_W;
        layers{end+1} = nnLinearLayer(W3, zeros(size(W3, 1)));

        W4 = dlt_layer.Operation_4_MatMul_W;
        layers{end+1} = nnLinearLayer(W4, zeros(size(W4, 1)));

        W5 = dlt_layer.Operation_5_MatMul_W;
        layers{end+1} = nnLinearLayer(W5, zeros(size(W5, 1)));

        W6 = dlt_layer.Operation_6_MatMul_W;
        layers{end+1} = nnLinearLayer(W6, zeros(size(W, 6)));

        W7 = dlt_layer.linear_7_MatMul_W;
        layers{end+1} = nnLinearLayer(W7, zeros(size(W, 7)));
    else
        % show warning
        if verbose
            warning("Skipping '%s'. Not implemented in cora yet!.", class(dlt_layer))
        end
        continue
    end

    currentSize = layers{end}.getOutputSize(currentSize);
end
layers = reshape(layers, [], 1); % 1 column

% instantiate neural network
obj = neuralNetwork(layers);

if ~isempty(inputSize)
    if isscalar(inputSize)
        inputSize = [inputSize, 1];
    end

    % set input size
    obj.setInputSize(inputSize, false);
    
    % sanity check (should not fail)
    x = reshape(zeros(inputSize), [], 1);
    obj.evaluate(x);
end

if verbose
    display(obj)
end

% ------------------------------ END OF CODE ------------------------------
