classdef nnCompositeLayer < nnLayer
% nnCompositeLayer - realize multiple parallel computation paths, e.g. res
%    connections
%
% Syntax:
%    obj = nnCompositeLayer(layers,aggregation)
%
% Inputs:
%    layers - cell array with the different parallel computation paths
%    aggregation - 'add' or 'concant'
%
% Outputs:
%    obj - generated object
%
% Example:
%   % A res connection.
%   obj = nnCompositeLayer({{nnLinearLayer(1,0), nnReLULayer}; {}}, 'add')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Lukas Koller
% Written:       27-June-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    layers
    aggregation
    outputSizes
    concatDim
end

methods
    % constructor
    function obj = nnCompositeLayer(layers, aggregation, varargin)
        obj@nnLayer(varargin{:})
        obj.layers = layers;
        obj.aggregation = aggregation;
        switch obj.aggregation
            case 'add'
                obj.aggregation = aggregation;
            case 'concat'
                obj.aggregation = aggregation;
            otherwise
                throw(CORAerror('CORA:wrongInputInConstructor', ...
                    'nnCompositeLayer.aggregation', ...
                    "Only supported value is 'add' and 'concat'!"));
        end

        % By default the concatenation is along the first dimension.
        obj.concatDim = 1;
    end

    function outputSize = getOutputSize(obj, inputSize)
        % Compute the output size of the first computation path.
        obj.outputSizes = cell(length(obj.layers),1);
        for i=1:length(obj.layers)
            layersi = obj.layers{i};
            obj.outputSizes{i} = inputSize;
            for j=1:length(layersi)
                obj.outputSizes{i} = ...
                    layersi{j}.computeSizes(obj.outputSizes{i});
            end
        end

        switch obj.aggregation
            case 'add'
                % Same output size as the individual outputs.
                outputSize = obj.outputSizes{1};
            case 'concat'
                outputSize = obj.outputSizes{1};
                outputSize_ = sum(cell2mat(obj.outputSizes),obj.concatDim);
                outputSize(obj.concatDim) = outputSize_(obj.concatDim);
            otherwise
                throw(CORAerror('CORA:wrongFieldValue', ...
                    'nnCompositeLayer.aggregation', ...
                    "Only supported value is 'add' and 'concat'!"));
        end
    end

    function [nin, nout] = getNumNeurons(obj)
        if isempty(obj.inputSize)
            nin = [];
            nout = [];
        else
            % We can only compute the number of neurons if the input
            % size was set.
            nin = prod(obj.inputSize);
            outputSize = getOutputSize(obj, obj.inputSize);
            nout = prod(outputSize);
        end
    end
end


% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})
    
    % numeric
    function r = evaluateNumeric(obj, input, options)
        obj.checkInputSize();

        % Initialize with neural element.
        r = [];
        imgSize = [];

        for i=1:length(obj.layers)
            % Compute result for the i-th computation path.
            layersi = obj.layers{i};
            ri = input;
            for j=1:length(layersi)        
                % Store input for backpropgation.
                if options.nn.train.backprop
                    layersi{j}.backprop.store.input = ri;
                end
                ri = layersi{j}.evaluateNumeric(ri, options);
            end
            % Obtain the output size of the current computation path.
            imgSizei = obj.outputSizes{i};
            
            if isempty(r)
                r = ri;
                imgSize = imgSizei;
            else
                switch obj.aggregation
                    case 'add'
                        % Add results.
                        r = r + ri;
                    case 'concat'
                        % Concatenate results.
                        [r,imgSize] = aux_concat(obj,r,imgSize,ri,imgSizei);
                    otherwise
                        throw(CORAerror('CORA:wrongFieldValue', ...
                            'nnCompositeLayer.aggregation', ...
                            "Only supported value is 'add' and 'concat'!"));
                end
            end
        end
    end

    % interval
    function r = evaluateInterval(obj, input, options)
        obj.checkInputSize()

        % Initialize with neural element.
        r = [];
        imgSize = [];

        for i=1:length(obj.layers)
            % Compute result for the i-th computation path.
            layersi = obj.layers{i};
            ri = input;
            for j=1:length(layersi)        
                % Store input for backpropgation.
                if options.nn.train.backprop
                    layersi{j}.backprop.store.input = ri;
                end
                ri = layersi{j}.evaluateInterval(ri, options);
            end
            % Obtain the output size of the current computation path.
            imgSizei = obj.outputSizes{i};
            
            if representsa(r,'emptySet')
                r = ri;
                imgSize = imgSizei;
            else
                switch obj.aggregation
                    case 'add'
                        % Add results.
                        r = r + ri;
                    case 'concat'
                        % Concatenate results.
                        [rl,~] = aux_concat(obj,r.inf,imgSize,ri.inf,imgSizei);
                        [ru,imgSize] = aux_concat(obj,r.sup,imgSize,ri.sup,imgSizei);
                        r = interval(rl,ru);
                    otherwise
                        throw(CORAerror('CORA:wrongFieldValue', ...
                            'nnCompositeLayer.aggregation', ...
                            "Only supported value is 'add' and 'concat'!"));
                end
            end
        end
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, options)
        obj.checkInputSize();

        % Retain input sensitivity.
        switch obj.aggregation
            case 'add'
                Sin = S;
            case 'concat'
                S_ = permute(S,[2 1 3]);
                Sin = aux_divideInput(obj,S_,obj.outputSizes);
                for i=1:length(Sin)
                    Sin{i} = permute(Sin{i},[2 1 3]);
                end
            otherwise
                throw(CORAerror('CORA:wrongFieldValue', ...
                    'nnCompositeLayer.aggregation', ...
                    "Only supported value is 'add' and 'concat'!"));
        end

        % Initialize with neural element.
        S = 0;

        for i=1:length(obj.layers)
            % Compute result for the i-th computation path.
            layersi = obj.layers{i};
            switch obj.aggregation
                case 'add'
                    Si = Sin;
                case 'concat'
                    Si = Sin{i};
                otherwise
                    throw(CORAerror('CORA:wrongFieldValue', ...
                        'nnCompositeLayer.aggregation', ...
                        "Only supported value is 'add' and 'concat'!"));
            end

            for j=length(layersi):-1:1
                Si = layersi{j}.evaluateSensitivity(Si,options);
            end
            % Aggregate results.
            S = S + Si;
        end

        if options.nn.store_sensitivity
            % Store the gradient (used for the sensitivity computation).
            obj.sensitivity = S;
        end
    end

    % zonotope batch
    function [rc, rG] = evaluateZonotopeBatch(obj, c, G, options)
        obj.checkInputSize();

        rc = [];
        rG = [];
        imgSize = [];

        for i=1:length(obj.layers)
            % Initialize output of the i-th computation path. 
            rci = c;
            rGi = G;
            % Compute result for the i-th computation path.
            layersi = obj.layers{i};
            for j=1:length(layersi)
                layersij = layersi{j};
                % Store input for backpropgation.
                if options.nn.train.backprop
                    layersij.backprop.store.inc = rci;
                    layersij.backprop.store.inG = rGi;
                end
                [rci,rGi] = layersij.evaluateZonotopeBatch(rci,rGi,options);
            end
            % Obtain the output size of the current computation path.
            imgSizei = obj.outputSizes{i};

            if isempty(rc)
                rc = rci;
                rG = rGi;
                imgSize = imgSizei;
            else
                switch obj.aggregation
                    case 'add'
                        % Add final results.
                        rc = rc + rci;
                        if size(rG,2) < size(rGi,2)
                            rGi(:,1:size(rG,2),:) = rGi(:,1:size(rG,2),:) + rG;
                            rG = rGi;
                        else
                            rG(:,1:size(rGi,2),:) = rG(:,1:size(rGi,2),:) + rGi;
                        end
                    case 'concat'
                        % Concatenate the centers.
                        [rc,~] = aux_concat(obj,rc,imgSize,rci,imgSizei);
                        % Concatenate generator matrices.
                        [rG,imgSize] = aux_concat(obj,rG,imgSize,rGi,imgSizei);
                    otherwise
                        throw(CORAerror('CORA:wrongFieldValue', ...
                            'nnCompositeLayer.aggregation', ...
                            "Only supported value is 'add' and 'concat'!"));
                end
            end
        end
    end

    function grad_in = backpropNumeric(obj, input, grad_out, options, updateWeights)
        obj.checkInputSize();

        switch obj.aggregation
            case 'add'
                %
            case 'concat'
                grad_out = aux_divideInput(obj,grad_out,obj.outputSizes);
            otherwise
                throw(CORAerror('CORA:wrongFieldValue', ...
                    'nnCompositeLayer.aggregation', ...
                    "Only supported value is 'add' and 'concat'!"));
        end

        % Initialize the gradient.
        grad_in = 0;
        for i=1:length(obj.layers)
            % Initialize output of the i-th computation path. 
            switch obj.aggregation
                case 'add'
                    grad_in_i = grad_out;
                case 'concat'
                    grad_in_i = grad_out{i};
                otherwise
                    throw(CORAerror('CORA:wrongFieldValue', ...
                        'nnCompositeLayer.aggregation', ...
                        "Only supported value is 'add' and 'concat'!"));
            end
            % Backpropagate the gradient through the i-th computation path.
            layersi = obj.layers{i};
            for j=length(layersi):-1:1
                % Retrieve stored input.
                inputsij = layersi{j}.backprop.store.input;
                % Compute the gradient.
                grad_in_i = layersi{j}.backpropNumeric(inputsij, ...
                    grad_in_i,options,updateWeights);
            end

            % Add final results.
            grad_in = grad_in + grad_in_i;
        end
    end

    % interval batch
    function [rgl, rgu] = backpropIntervalBatch(obj, l, u, gl, gu, options, updateWeights)
        switch obj.aggregation
            case 'add'
                %
            case 'concat'
                gl = aux_divideInput(obj,gl,obj.outputSizes);
                gu = aux_divideInput(obj,gu,obj.outputSizes);
            otherwise
                throw(CORAerror('CORA:wrongFieldValue', ...
                    'nnCompositeLayer.aggregation', ...
                    "Only supported value is 'add' and 'concat'!"));
        end

        rgl = 0;
        rgu = 0;
        for i=1:length(obj.layers)
            % Initialize output of the i-th computation path. 
            switch obj.aggregation
                case 'add'
                    rgli = gl;
                    rgui = gu;
                case 'concat'
                    rgli = gl{i};
                    rgui = gu{i};
                otherwise
                    throw(CORAerror('CORA:wrongFieldValue', ...
                        'nnCompositeLayer.aggregation', ...
                        "Only supported value is 'add' and 'concat'!"));
            end
            % Compute result for the i-th computation path.
            layersi = obj.layers{i};
            idxLayeri = flip(1:length(layersi));
            for j=idxLayeri
                layersij = layersi{j};
                % Retrieve stored input
                input = layersij.backprop.store.input;
                l = input.inf;
                u = input.sup;
                % Compute the gradient.
                [rgli,rgui] = layersij.backpropIntervalBatch(l,u,rgli,rgui, ...
                    options,updateWeights);
            end

            % Add final results.
            rgl = rgl + rgli;
            rgu = rgu + rgui;
        end
    end

    % zonotope batch
    function [rgc, rgG] = backpropZonotopeBatch(obj, c, G, ...
            gc, gG, options, updateWeights)
        switch obj.aggregation
            case 'add'
                %
            case 'concat'
                gc = aux_divideInput(obj,gc,obj.outputSizes);
                gG = aux_divideInput(obj,gG,obj.outputSizes);
            otherwise
                throw(CORAerror('CORA:wrongFieldValue', ...
                    'nnCompositeLayer.aggregation', ...
                    "Only supported value is 'add' and 'concat'!"));
        end

        rgc = 0;
        rgG = 0;
        for i=1:length(obj.layers)
            % Initialize output of the i-th computation path. 
            switch obj.aggregation
                case 'add'
                    rgci = gc;
                    rgGi = gG;
                case 'concat'
                    rgci = gc{i};
                    rgGi = gG{i};
                otherwise
                    throw(CORAerror('CORA:wrongFieldValue', ...
                        'nnCompositeLayer.aggregation', ...
                        "Only supported value is 'add' and 'concat'!"));
            end
            % Compute result for the i-th computation path.
            layersi = obj.layers{i};
            idxLayeri = flip(1:length(layersi));
            for j=idxLayeri
                layersij = layersi{j};
                % Retrieve stored input
                c = layersij.backprop.store.inc;
                G = layersij.backprop.store.inG;
                % Compute the gradient.
                [rgci,rgGi] = layersij.backpropZonotopeBatch(c,G, ...
                    rgci,rgGi,options,updateWeights);
            end

            % Add final results.
            rgc = rgc + rgci;
            if size(rgGi,2) < size(rgG,2)
                rgG(:,1:size(rgGi,2),:) = rgG(:,1:size(rgGi,2),:) + rgGi;
            else
                rgGi(:,1:size(rgG,2),:) = rgGi(:,1:size(rgG,2),:) + rgG;
                rgG = rgGi;
            end
        end
    end
end

% Auxiliary functions -----------------------------------------------------

methods

    function [x,imgSize] = aux_concat(obj,x1,imgSize1,x2,imgSize2)

        % Check if we are concatenating matrices, e.g., generator matrix,
        % sensitivity.
        isMatrixConcat = ndims(x1) > 2;
        if isMatrixConcat
            % Reshape the results for easier concatenation; we move the
            % extra dimension to the batch.
            [n1,q1,bSz] = size(x1);
            [n2,q2,~] = size(x2);
            % Pad the generator matrices with zeros.
            if q1 < q2
                x1 = cat(2,x1,zeros([n1 q2-q1 bSz],'like',x1));
            elseif q2 < q1
                x2 = cat(2,x2,zeros([n2 q1-q2 bSz],'like',x2));
            end
            % Update the number of generators.
            q = max(q1,q2);
            % Reshape the matrices.
            x1 = x1(:,:);
            x2 = x2(:,:);
        end
        % Obtain the batch size.
        [~,bSz_] = size(x1);

        % Reshape the inputs.
        x1_ = reshape(x1,[imgSize1 bSz_]);
        x2_ = reshape(x2,[imgSize2 bSz_]);
        % Concatenate the images.
        x_ = cat(obj.concatDim,x1_,x2_);
        % Reshape the result.
        x = reshape(x_,[],bSz_);
        % Compute concatenated image size.
        imgSize = imgSize1;
        imgSize(obj.concatDim) = imgSize(obj.concatDim) ...
            + imgSize2(obj.concatDim);

        if isMatrixConcat
            x = reshape(x,[],q,bSz);
        end
    end

    function xis = aux_divideInput(obj,x,imgSizes)
        % Initialize the indices.
        idx = [];
        imgSize = [];
        for i=1:length(imgSizes)
            % Obtain the i-th image size.
            imgSizei = imgSizes{i};
            % Create indices.
            idxi = i*ones(imgSizei);
            % Concatenate the indices.
            if isempty(idx)
                idx = idxi;
                imgSize = imgSizei;
            else
                [idx,imgSize] = aux_concat(obj,idx,imgSize,idxi,imgSizei);
            end
        end

        % Check if we are concatenating matrices, e.g., generator matrix,
        % sensitivity.
        isMatrixConcat = ndims(x) > 2;
        if isMatrixConcat
            % Reshape the results for easier concatenation; we move the
            % extra dimension to the batch.
            [~,q,bSz] = size(x);
            x = x(:,:);
        end

        % Initialize the results.
        xis = cell(length(imgSizes),1);
        for i=1:length(imgSizes)
            % Extract the elements.
            xis{i} = x(idx == i,:);
            
            if isMatrixConcat
                xis{i} = reshape(xis{i},[],q,bSz);
            end
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
