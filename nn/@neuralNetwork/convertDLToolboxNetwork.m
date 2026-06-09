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

% Authors:       Tobias Ladner, Lukas Koller, Benedikt Kellner
% Written:       30-March-2022
% Last update:   05-June-2022 (LK, Conv, Pool)
%                17-January-2023 (TL, Reshape)
%                23-November-2023 (TL, bug fix with scalar element-wise operation)
%                25-July-2023 (TL, nnElementwiseAffineLayer)
%                31-July-2023 (LK, nnSoftmaxLayer)
%                14-March-2026 (BK, ScalingLayer, CustomOutputLayer)
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
        
        try
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
        catch ME
            if iscell(dlt_layer)
                 name = "Composite";
                 cls = "cell";
            else
                 name = dlt_layer.Name;
                 cls = class(dlt_layer);
            end
            fprintf(2, "Error converting layer %d (%s, class: %s):\n%s\n", ...
                i, name, cls, ME.message);
            rethrow(ME);
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
    
    elseif isa(dlt_layer, 'nnet.onnx.layer.ElementwiseAffineLayer') || ...
            isa(dlt_layer, 'nnet.cnn.layer.ScalingLayer')
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

         % sanity check
        if length(size(W)) == 4
             assert(size(W, 3) == dlt_layer.NumChannels, ...
            'CORA:convertDLToolboxNetwork:Conv2DDimensions', ...
            'Expected 3rd dimension of Weights to match NumChannels (%d), but got %d.', ...
            dlt_layer.NumChannels, size(W, 3));
             assert(size(W, 4) == dlt_layer.NumFilters, ...
            'CORA:convertDLToolboxNetwork:Conv2DDimensions', ...
            'Expected 4th dimension of Weights to match NumFilters (%d), but got %d.', ...
            dlt_layer.NumFilters, size(W, 4));
        end

    elseif isa(dlt_layer, 'nnet.cnn.layer.TransposedConvolution2DLayer')
        W = double(dlt_layer.Weights);
        b = double(reshape(dlt_layer.Bias, [], 1));
        cs = dlt_layer.CroppingSize;   % DLT: [top bottom left right]
    
        if numel(cs) == 4
            % Convert DLT [top bottom left right] -> CORA [left top right bottom]
            cropping = [cs(3), cs(1), cs(4), cs(2)];
        else
            % handle scalar / 1x2 / 2x2 cases as needed
            cs = cs(:)';
            cropping = [cs(2), cs(1), cs(2), cs(1)];
        end
    
        stride = dlt_layer.Stride;
        dilation = [1, 1];  % DLT transposed conv has no dilation
    
        layers{end+1} = nnConvTranspose2DLayer(W, b, cropping, stride, dilation, dlt_layer.Name);

        % sanity check
        if length(size(W)) == 4
             assert(size(W, 3) == dlt_layer.NumFilters, ...
            'CORA:convertDLToolboxNetwork:TransposedConvDimensions', ...
            'Expected 3rd dimension of Weights to match NumFilters (%d), but got %d.', ...
            dlt_layer.NumFilters, size(W, 3));
             assert(size(W, 4) == dlt_layer.NumChannels, ...
            'CORA:convertDLToolboxNetwork:TransposedConvDimensions', ...
            'Expected 4th dimension of Weights to match NumChannels (%d), but got %d.', ...
            dlt_layer.NumChannels, size(W, 4));
        end

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

    elseif contains(class(dlt_layer), 'Upsample')
        % upsample
        
        % Check for interpolation mode if available (e.g. in ONNX layers)
        if isprop(dlt_layer, 'Mode')
            mode = dlt_layer.Mode;
            if (isa(mode, 'string') || isa(mode, 'char')) && ~strcmpi(mode, 'nearest')
                throw(CORAerror('CORA:notDefined', ...
                      'convertDLToolboxNetwork/UpsampleMode: Only "nearest" neighbor upsampling is supported. Found mode: %s', mode));
            end
        end
        
        % determine scale factor
        scale = [];
         
         % Check properties starting with onnx__
         if isempty(scale)
             props = properties(dlt_layer);
             for k=1:length(props)
                 if contains(props{k}, 'onnx__')
                     val = dlt_layer.(props{k});
                     try
                         if isa(val, 'dlarray'); val = extractdata(val); end
                         val = double(val);
                     catch
                          props = properties(dlt_layer);
                          throw(CORAerror('CORA:wrongFieldValue', ...
                             'convertDLToolboxNetwork/UpsampleScale: Could not detect Scale/Scales for Upsample layer. Properties: %s', strjoin(props, ', '))); 
                     end
                     if numel(val)==4 && all(val(1:2)==1) && all(val(3:4)>1)
                         scale = val; break;
                     end
                 end
             end
         end

        if length(scale) == 4
            scale = scale(3:4); 
        end
        
        if isscalar(scale)
             scale = [scale scale];
        end
        
        currentSize = [currentSize(:)', ones(1, max(0, 3-length(currentSize)))];
        H = currentSize(1);
        W = currentSize(2);
        C = currentSize(3);
        
        % Nearest Neighbor reindexing
        % generate input indices
        rows = repelem((1:H)', scale(1));
        cols = repelem(1:W, scale(2));
        
        % mapping: output(i,j) -> input(rows(i), cols(j))
        % linear index in input: (cols-1)*H + rows
        linIdx = rows + (cols-1)*H;
        
        % handle channels
        offsets = reshape((0:C-1)*H*W, 1, 1, C);
        idx_out = linIdx + offsets;
        
        layers{end+1} = nnReshapeLayer(idx_out, dlt_layer.Name);

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
            isa(dlt_layer, 'nnet.cnn.layer.ClassificationOutputLayer') || ...
            isa(dlt_layer, 'nnet.onnx.layer.CustomOutputLayer')
        % ignore
        inputSize = [];
        return

    % Custom Layers -----------------------------------------------
    elseif strcmp(dlt_layer.Name, 'Gemm_To_ReshapeLayer1000')
        % TODO: safe as a CORA linear layer, followed by a reshpae layer
        % Read the corresponding weights from the dl toolbox layer
        W = double(dlt_layer.gen_l1_0_weight');   % now 512 x 5
        b = double(dlt_layer.gen_l1_0_bias);      % 512 x 1
        layers{end+1} = nnLinearLayer(W, b, dlt_layer.Name);

        % Get all field names of the 'Vars' structure
        varNames = fieldnames(dlt_layer.Vars);
        
        % Find the field name that contains the specific string
        matchingFieldIndices = contains(varNames, 'onnx__Reshape_');
        
        % Check if a match was found and get the name (assuming only one match)
        if any(matchingFieldIndices)
            fieldName = varNames{matchingFieldIndices};
            
            % Use dynamic field access
            onnxShape = double(extractdata(dlt_layer.Vars.(fieldName)));
        else
            % Handle case where no field matches the pattern
            CORAwarning('CORA:nn','No field containing "onnx__Reshape_" found in dlt_layer.Vars.');
            onnxShape = []; % or some other default/error value
        end

        % Convert to DLT/MATLAB ordering and drop batch dim:
        shapeDlt = flip(onnxShape);          % [2;2;128;1]
        H = shapeDlt(1);
        W = shapeDlt(2);
        C = shapeDlt(3);
        
        % Build index tensor with desired output shape [H W C]
        outSize = [H, W, C];                 % [2 2 128]
        idx_out = reshape(1:prod(outSize), outSize);
        
        % CORA reshape layer equivalent to your ONNX Reshape/Flatten
        layers{end+1} = nnReshapeLayer(idx_out, dlt_layer.Name);
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
        layers{end+1} = nnLinearLayer(params.Learnables.linear_7_MatMul_W, params.Nonlearnables.linear_7_Add_B,dlt_layer.Name)
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
