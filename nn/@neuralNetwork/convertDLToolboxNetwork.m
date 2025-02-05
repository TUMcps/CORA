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
%                25-July-2023 (TL, nnElementwiseAffineLayer)
%                31-July-2023 (LK, nnSoftmaxLayer)
% Last revision: 17-August-2022

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 2
    verbose = false;
end

n = length(dltoolbox_layers);

layers = {};
inputSize = [];
currentSize = [];

if verbose
    disp("Converting Deep Learning Toolbox Model to neuralNetwork...")
end

for i = 1:n
    dlt_layer = dltoolbox_layers{i};
    if verbose
        fprintf("#%d: %s\n", i, class(dlt_layer))
    end

    if iscell(dlt_layer)
        % We need to construct a composite layer.
        layers_ = {};
        currentSize_ = currentSize;
        for j=1:length(dlt_layer)
            layersj_ = {};
            for k=1:length(dlt_layer{j})
                [layersj_,~,currentSize_] = aux_convertLayer(layersj_, ...
                    dlt_layer{j}(k),currentSize_,verbose);
            end
            layers_{j} = layersj_;
        end
        % Obtain aggregation layer
        i = i+1;
        aggr_dlt_layer = dltoolbox_layers{i};
        if isa(aggr_dlt_layer,'nnet.cnn.layer.AdditionLayer')
            aggregation = 'add';
        else
            % Error
            aggregation = [];
        end
        compLayer = nnCompositeLayer(layers_,aggregation);
        layers{end+1} = compLayer;

        currentSize = layers{end}.getOutputSize(currentSize);
    else
        % Just append a regular layer.
        [layers,inputSize_,currentSize] = aux_convertLayer(layers, ...
            dlt_layer,currentSize,verbose);
        if isempty(inputSize)
            inputSize = inputSize_;
        end
    end
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


% Auxiliary functions -----------------------------------------------------

function [layers,inputSize,currentSize] = aux_convertLayer(layers,dlt_layer,currentSize,verbose)
    % handle different types of layers
    if isa(dlt_layer, 'nnet.cnn.layer.ImageInputLayer') || ...
            isa(dlt_layer, 'nnet.cnn.layer.FeatureInputLayer')
        inputSize = dlt_layer.InputSize;
        if length(inputSize) == 1
            inputSize = [inputSize, 1];
        elseif length(inputSize) == 3
            % channel dimension should be last: [h,w,c]
            % inputSize = sort(inputSize, 'descend');
        end

        if strcmp(dlt_layer.Normalization,'zscore')
            mu = dlt_layer.Mean;
            sigma = dlt_layer.StandardDeviation;
            layers{end+1} = nnElementwiseAffineLayer(1/sigma, ...
                -mu/sigma,dlt_layer.Name);
        end
        currentSize = inputSize;
        return

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
                else
                    s = reshape(repmat(s,currentSize./size(s)),currentSize);
                end
            end
            if ~isscalar(o)
                % fix if all values are equal
                if ~isempty(o) && all(o(1) == o,'all')
                    o = o(1);
                % try to fix offset vector
                else
                    o = reshape(repmat(o,currentSize./size(o)),currentSize);
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
        % average pooling 
        poolSize = dlt_layer.PoolSize;
        padding = dlt_layer.PaddingSize;
        stride = dlt_layer.Stride;
        dilation = [1, 1];

        layers{end+1} = nnAvgPool2DLayer(poolSize, padding, stride, dilation, dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.Convolution2DLayer')
        % convolutional 2D
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

        % activation functions ---
    elseif isa(dlt_layer, 'nnet.cnn.layer.ReLULayer')
        % relu
        layers{end+1} = nnReLULayer(dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.LeakyReLULayer')
        % leaky relu
        alpha = double(dlt_layer.Scale);
        layers{end+1} = nnLeakyReLULayer(alpha, dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.TanhLayer')
        % tanh
        layers{end+1} = nnTanhLayer(dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.SigmoidLayer')
        % sigmoid
        layers{end+1} = nnSigmoidLayer(dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.SoftmaxLayer')
        % softmax
        layers{end+1} = nnSoftmaxLayer(dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.onnx.layer.IdentityLayer')
        % ignore
        inputSize = [];
        return

    elseif isa(dlt_layer, 'nnet.cnn.layer.RegressionOutputLayer') || ...
            isa(dlt_layer, 'nnet.cnn.layer.ClassificationOutputLayer')
        % ignore
        inputSize = [];
        return

        % Custom Layers -----------------------------------------------

    elseif contains(lower(class(dlt_layer)), 'flatten') || ...
            contains(lower(class(dlt_layer)), 'reshape')
        % flatten
        
        idx = dlarray(1:prod(currentSize));
        idx = reshape(idx, currentSize);
        idx_out = dlt_layer.predict(idx);
        idx_out = double(extractdata(idx_out));
        layers{end+1} = nnReshapeLayer(idx_out, dlt_layer.Name);

    elseif strcmp(dlt_layer.Name, 'MatMul_To_ReluLayer1003')
        % for VNN Comp (test -- test_nano.onnx)
        params = dlt_layer.ONNXParams;
        layers{end+1} = nnLinearLayer(params.Learnables.Ma_MatMulcst,0, ...
            dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
    elseif strcmp(dlt_layer.Name, 'MatMul_To_AddLayer1003')
        % for VNN Comp (test)
        params = dlt_layer.ONNXParams;
        if ~isfield(params.Learnables,'W2')
            % (test --- test_tiny.onnx)
            layers{end+1} = nnLinearLayer(params.Learnables.W0,0,dlt_layer.Name);
            layers{end+1} = nnReLULayer(dlt_layer.Name);
            layers{end+1} = nnLinearLayer(params.Learnables.W1,0,dlt_layer.Name);
        else
            % (test --- test_small.onnx)
            layers{end+1} = nnLinearLayer(params.Learnables.W0',[1.5; 1.5],dlt_layer.Name);
            layers{end+1} = nnReLULayer(dlt_layer.Name);
            layers{end+1} = nnLinearLayer(params.Learnables.W1,[2.5; 2.5],dlt_layer.Name);
            layers{end+1} = nnReLULayer(dlt_layer.Name);
            layers{end+1} = nnLinearLayer(params.Learnables.W2',3.5,dlt_layer.Name);
        end
    elseif strcmp(dlt_layer.Name, 'MatMul_To_AddLayer1019') || ...
            strcmp(dlt_layer.Name, 'Mul_To_AddLayer1021')
        % for VNN Comp (cora - mnist, svhn, cifar10)
        % requires hard-coding ...
        params = dlt_layer.ONNXParams;
        layers{end+1} = nnLinearLayer(params.Learnables.fc_1_copy_MatMul_W, params.Nonlearnables.fc_1_copy_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.fc_2_copy_MatMul_W, params.Nonlearnables.fc_2_copy_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.fc_3_copy_MatMul_W, params.Nonlearnables.fc_3_copy_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.fc_4_copy_MatMul_W, params.Nonlearnables.fc_4_copy_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.fc_5_copy_MatMul_W, params.Nonlearnables.fc_5_copy_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.fc_6_copy_MatMul_W, params.Nonlearnables.fc_6_copy_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.fc_7_copy_MatMul_W, params.Nonlearnables.fc_7_copy_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.fc_8_copy_MatMul_W, params.Nonlearnables.fc_8_copy_Add_B,dlt_layer.Name);

    elseif strcmp(dlt_layer.Name, 'Sub_To_AddLayer1018')
        % for VNN Comp (test_sat.onnx)
        % requires hard-coding ...
        params = dlt_layer.ONNXParams;
        layers{end+1} = nnLinearLayer(params.Learnables.Operation_1_MatMul_W, params.Nonlearnables.Operation_1_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.Operation_2_MatMul_W, params.Nonlearnables.Operation_2_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.Operation_3_MatMul_W, params.Nonlearnables.Operation_3_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.Operation_4_MatMul_W, params.Nonlearnables.Operation_4_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.Operation_5_MatMul_W, params.Nonlearnables.Operation_5_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.Operation_6_MatMul_W, params.Nonlearnables.Operation_6_Add_B,dlt_layer.Name);
        layers{end+1} = nnReLULayer(dlt_layer.Name);
        layers{end+1} = nnLinearLayer(params.Learnables.linear_7_MatMul_W, params.Nonlearnables.linear_7_Add_B,dlt_layer.Name);
    else
        % unknown layer, show warning
        if verbose
            CORAwarning('CORA:nn',"Skipping '%s'. Not implemented in cora yet!.", class(dlt_layer))
        end
        inputSize = [];
        return
    end
    currentSize = layers{end}.getOutputSize(currentSize);
    inputSize = [];
end

end

% ------------------------------ END OF CODE ------------------------------
