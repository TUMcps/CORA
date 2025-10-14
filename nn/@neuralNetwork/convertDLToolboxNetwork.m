function obj = convertDLToolboxNetwork(dlt_layers, verbose)
% convertDLToolboxNetwork - converts a network from the Deep Learning
%    Toolbox to a CORA neuralNetwork for verification
%
% Syntax:
%    res = neuralNetwork.convertDLToolboxNetwork(dlt_layers)
%    res = neuralNetwork.convertDLToolboxNetwork(dlt_layers,verbose)
%
% Inputs:
%    dlt_layers - layer array (e.g. dltoolbox_nn.Layers)
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

if verbose
    disp("Converting Deep Learning Toolbox Model to neuralNetwork...")
end

% Initialize input size and current size.
inputSize = [];
currentSize = [];

% Convert the layers in the nested cell array.
[layers,inputSize,~] = aux_convertLayers(dlt_layers, ...
    inputSize,currentSize,verbose);

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

function [layers,inputSize,currentSize] = aux_convertLayers( ...
        dlt_layers,inputSize,currentSize,verbose)
    % Initialize the layers cell array.
    layers = {};

    % Needed to separate different input for reshape layers that are
    % followed by a composite layer.
    nextInputIdx = {};

    % Recursively convert the layers in the nested cell array.
    for i=1:length(dlt_layers)
        % Obtain the i-th layer.
        dlt_layer = dlt_layers{i};
        
        if iscell(dlt_layer)
            if verbose
                fprintf("#%d: Composite\n", i)
            end
            % We need to construct a composite layer.
            layersi = {};
            % Iterate the computation paths.
            for j=1:length(dlt_layer)
                if ~isempty(nextInputIdx)
                    % Construct the reshape layers that reshapes the input
                    % of the j-th computation path.
                    reshapeIdx = reshape(nextInputIdx{j},1,[]);
                    preshapeLayer = nnReshapeLayer( ...
                        reshapeIdx,dlt_layers{i-1}.Name);
                    currentSize = size(nextInputIdx{j});
                end
                % Convert the j-th computation path.
                [layersij,~,~] = aux_convertLayers( ...
                    dlt_layer{j},inputSize,currentSize,verbose);
                if ~isempty(nextInputIdx)
                    % Prepend the reshape layer.
                    layersij = [{preshapeLayer}; layersij];
                end
                % Append the converted layers.
                layersi{j} = layersij;
            end
            % Increment the index to obtain aggregation layer.
            i = i+1;
            % Obtain the aggregation layer.
            aggr_dlt_layer = dlt_layers{i};
            % Check the type of aggregation layer.
            if isa(aggr_dlt_layer,'nnet.cnn.layer.AdditionLayer')
                aggregation = 'add';
            elseif isa(aggr_dlt_layer,'nnet.cnn.layer.ConcatenationLayer')
                aggregation = 'concat';
            else
                % Aggregation type is not supported.
                aggregation = [];
            end
            % Instantiate the composite layer.
            compLayer = nnCompositeLayer(layersi,aggregation);
            % Append the composite layer.
            layers{end+1} = compLayer;
            % Update the output size.
            currentSize = layers{end}.getOutputSize(currentSize);
            % Reset the input indices.
            nextInputIdx = {};
        else
            if verbose
                fprintf("#%d: %s\n", i, class(dlt_layer))
            end
            % Just append a regular layer.
            [layers,inputSize_,currentSize,nextInputIdx] = ...
                aux_convertLayer(layers,dlt_layer,currentSize,verbose);
            if isempty(inputSize)
                inputSize = inputSize_;
            end
        end
    end
    layers = reshape(layers, [], 1); % 1 column
end

function [layers,inputSize,currentSize,nextInputIdx] = ...
        aux_convertLayer(layers,dlt_layer,currentSize,verbose)
    % By default all inputs are given to the next layer; for reshape layers
    % with multiple input we can encode the path of the inputs.
    nextInputIdx = {};

    % Initialize the input size.
    inputSize = [];

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
                -mu./sigma,dlt_layer.Name);
        end
        currentSize = inputSize;
        return
    elseif isa(dlt_layer,'nnet.cnn.layer.SequenceInputLayer')
        inputSize = [dlt_layer.InputSize dlt_layer.MinLength];
        currentSize = inputSize;
        return

        % Normal Layers ---------------------------------------------------

        % linear layers ---

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

        % convolutional ---

    elseif isa(dlt_layer, 'nnet.cnn.layer.Convolution2DLayer')
        % convolutional 2D
        W = double(dlt_layer.Weights);
        b = double(reshape(dlt_layer.Bias, [], 1));
        padding = dlt_layer.PaddingSize;
        stride = dlt_layer.Stride;
        dilation = dlt_layer.DilationFactor;
        layers{end+1} = nnConv2DLayer(W, b, padding, stride, dilation, dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.BatchNormalizationLayer')
        % batch normalization
        mean = dlt_layer.TrainedMean;
        var = dlt_layer.TrainedVariance;
        epsilon = dlt_layer.Epsilon;
        scale = dlt_layer.Scale;
        bias = dlt_layer.Offset;
        layers{end+1} = nnBatchNormLayer(scale, bias, var, mean, dlt_layer.Name, epsilon);
        
    elseif isa(dlt_layer, 'nnet.cnn.layer.AveragePooling2DLayer')
        % average pooling 
        poolSize = dlt_layer.PoolSize;
        padding = dlt_layer.PaddingSize;
        stride = dlt_layer.Stride;
        dilation = [1, 1];

        layers{end+1} = nnAvgPool2DLayer(poolSize, padding, stride, dilation, dlt_layer.Name);

    elseif isa(dlt_layer, 'nnet.cnn.layer.MaxPooling2DLayer')
        % max pooling
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
            contains(lower(class(dlt_layer)), 'reshape') || ...
            contains(lower(class(dlt_layer)), 'slice')
        % flatten
        
        idx = dlarray(1:prod(currentSize));
        idx = reshape(idx, currentSize);

        idx_out = cell(1,dlt_layer.NumOutputs);
        [idx_out{:}] = dlt_layer.predict(idx);
        idx_out = cellfun(@(idx) double(extractdata(idx)),idx_out, ...
            'UniformOutput',false);

        if length(idx_out) > 1
            % There are multiple successor layer. We prepend the reshape
            % layer to each computation path.
            nextInputIdx = idx_out;
            return;
        else
            % There is only a single successor layer.
            layers{end+1} = nnReshapeLayer(idx_out{1}, dlt_layer.Name);
        end

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
    % Update the current size.
    currentSize = layers{end}.getOutputSize(currentSize);
end

end

% ------------------------------ END OF CODE ------------------------------
