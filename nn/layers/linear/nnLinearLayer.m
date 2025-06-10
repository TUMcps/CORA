classdef nnLinearLayer < nnLayer
% nnLinearLayer - class for linear layers
%
% Syntax:
%    obj = nnLinearLayer(W, b)
%    obj = nnLinearLayer(W, b, name)
%
% Inputs:
%    W - weight matrix
%    b - bias column vector
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner, Lukas Koller
% Written:       28-March-2022
% Last update:   23-November-2022 (polish)
%                14-December-2022 (variable input tests, inputArgsCheck)
%                03-May-2023 (LK, added backprop for polyZonotope)
%                25-May-2023 (LK, modified sampling of gradient for 'extreme')
%                25-July-2023 (LK, sampling of gradient with cartProd)
%                02-August-2023 (LK, added zonotope batch-eval & -backprop)
%                19-August-2023 (LK, zonotope batch-eval: memory optimizations for GPU training)
%                22-January-2024 (LK, functions for IBP-based training)
% Last revision: 10-August-2022 (renamed)

% ------------------------------ BEGIN CODE -------------------------------

properties
    W                       % weight matrix
    b                       % bias
    d = []                  % approx error (additive)
end

properties (Constant)
    is_refinable = false    % whether the layer is refineable
end

methods
    % constructor
    function obj = nnLinearLayer(W, varargin)
        % parse input
        [b, name] = setDefaultValues({0, []}, varargin);
        inputArgsCheck({ ...
            {W, 'att', {'numeric', 'gpuArray'}}; ...
            {b, 'att', {'numeric', 'gpuArray'}}; ...
        })

        % check dimensions
        if length(b) == 1
            b = b * ones(size(W, 1), 1);
        end
        if ~all(size(b, 1) == size(W, 1))
           throw(CORAerror('CORA:wrongInputInConstructor', ...
               'The dimensions of W and b should match.'));
        end
        if size(b, 2) ~= 1
           throw(CORAerror('CORA:wrongInputInConstructor', ...
               "Second input 'b' should be a column vector."));
        end

        % call super class constructor
        obj@nnLayer(name)

        obj.W = double(W);
        obj.b = double(b);
    end
end

% evaluate ----------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})
    
    % numeric
    function r = evaluateNumeric(obj, input, options)
        % linear transformation
        if representsa_(input,'emptySet',eps)
            r = obj.b .* ones(1,size(input,2));
        else
            r = obj.W * input + obj.b;
        end

        % add approx error
        if ~representsa_(obj.d,'emptySet',eps)
            samples = obj.d.randPoint(size(r,2));
            r = r + samples;
        end
    end

    % interval 
    function bounds = evaluateInterval(obj, bounds, options)
        % IBP (see Gowal et al. 2019)
        % Compute center and radius.
        mu = (bounds.sup + bounds.inf)/2;
        r = (bounds.sup - bounds.inf)/2;
        % Compute linear relaxation.
        mu = pagemtimes(obj.W,mu) + obj.b;
        r = pagemtimes(abs(obj.W),r);
        % Convert center and radius back to lower and upper bound.
        bounds = interval(mu - r,mu + r);

        % add approx error
        if ~representsa_(obj.d,'emptySet',eps)
            bounds = bounds + obj.d;
        end
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, options)
        % S = S * obj.W;
        % use pagemtimes to compute sensitivity simultaneously for an
        % entire batch.
        S = pagemtimes(S,obj.W);
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
        c = obj.W * c + obj.b;
        G = obj.W * G;
        GI = obj.W * GI;

        if ~representsa_(obj.d,'emptySet',eps)
            c = c + center(obj.d);
            GI = [GI, diag(rad(obj.d))];
        end
    end
    
    % zonotope batch (for training)
    function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
        [n,~,batchSize] = size(G);
        if options.nn.interval_center
            cl = reshape(c(:,1,:),[n batchSize]);
            cu = reshape(c(:,2,:),[n batchSize]);
            c = obj.evaluateInterval(interval(cl,cu));
            c = permute(cat(3,c.inf,c.sup),[1 3 2]);
        else
            c = obj.W*c + obj.b;
        end
        G = pagemtimes(obj.W,G);
    end

    % taylm
    function r = evaluateTaylm(obj, input, options)
        r = obj.W * input + obj.b;
    end

    % conZonotope
    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options)
        c = obj.W * c + obj.b;
        G = obj.W * G;
    end

    % backprop ------------------------------------------------------------

    function grad_in = backpropNumeric(obj, input, grad_out, options)
        % update weights and bias
        obj.updateGrad('W', grad_out * input', options);
        obj.updateGrad('b', sum(grad_out, 2), options);
        % backprop gradient
        grad_in = obj.W' * grad_out;
    end

    function [gl, gu] = backpropIntervalBatch(obj, l, u, gl, gu, options)
        % see Gowal et al. 2019
        mu = (u + l)/2;
        r = (u - l)/2;

        % update weights and bias
        obj.updateGrad('W', (gu + gl)*mu' + (gu - gl)*r'.*sign(obj.W), options);
        obj.updateGrad('b', sum(gu + gl, 2), options);

        % backprop gradient
        dmu = obj.W' * (gu + gl)/2;
        dr = abs(obj.W') * (gu - gl)/2;
        gl = dmu - dr;
        gu = dmu + dr;
    end

    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options)
        [n,numGen,batchSize] = size(G);
        % obtain indices of active generators
        genIds = obj.backprop.store.genIds;

        if strcmp(options.nn.train.zonotope_weight_update,'center')
            % use the center to update the weights and biases
            weightsUpdate = gc*c';
            biasUpdate = sum(gc,2);
        elseif strcmp(options.nn.train.zonotope_weight_update,'sample')
            % sample random point factors
            beta = 2*rand(numGen,1,batchSize,'like',c) - 1;
            % compute gradient samples
            grads = gc + reshape(pagemtimes(gG,beta),size(c));
            % compute input samples
            inputs = inc + reshape(pagemtimes(G,beta),size(c));
            % Compute weights and bias update
            weightsUpdate = grads*inputs';
            biasUpdate = sum(grads,2);
        elseif strcmp(options.nn.train.zonotope_weight_update,'extreme')
            numSamples = 1;
            % sample a point that has only factors {-1,1}
            beta = randi([-1,1],numGen,numSamples,batchSize,'like',c);
            % compute gradient samples
            grads = permute(repmat(gc,1,1,numSamples),[1 3 2]) + ...
                pagemtimes(gG,beta);
            % compute input samples
            inputs = permute(repmat(c,1,1,numSamples),[1 3 2]) + ...
                pagemtimes(G,beta);
            % Compute weights and bias update
            weightsUpdate = squeeze(mean(pagemtimes(grads,'none',...
                inputs,'transpose'),3));
            biasUpdate = squeeze(sum(mean(grads,2),3));
        elseif strcmp(options.nn.train.zonotope_weight_update,'outer_product')
            % compute outer product of gradient and input zonotope
            % (1) outer product between centers
            centerTerm = gc*c';
            % (2) outer product between generator matrices
            gensTerm = 1/3*sum(pagemtimes(gG(:,genIds,:),'none', ...
                G(:,genIds,:),'transpose'),3);
            % Compute weights and bias update
            weightsUpdate = centerTerm + gensTerm; 
            biasUpdate = sum(gc,2);
        elseif strcmp(options.nn.train.zonotope_weight_update,'sum')
            % compute outer product of gradient and input zonotope
            if options.nn.interval_center
                % (1) outer product between centers
                cl = reshape(c(:,1,:),[n batchSize]);
                cu = reshape(c(:,2,:),[n batchSize]);
                gl = reshape(gc(:,1,:),[size(gc,1) batchSize]);
                gu = reshape(gc(:,2,:),[size(gc,1) batchSize]);
                [gl, gu] = backpropIntervalBatch(obj, cl, cu, gl, gu, options);
                gc = permute(cat(3,gl,gu),[1 3 2]);
                % (2) outer product between generator matrices
                gensTerm = sum(pagemtimes(gG(:,genIds,:),'none', ...
                    G(:,genIds,:),'transpose'),3);
                % Compute weights and bias update
                weightsUpdate = gensTerm;
                biasUpdate = 0;
            else
                % (1) outer product between centers
                centerTerm = gc*c';
                % (2) outer product between generator matrices
                gensTerm = sum(pagemtimes(gG(:,genIds,:),'none', ...
                    G(:,genIds,:),'transpose'),3);
                % Compute weights and bias update
                weightsUpdate = centerTerm + gensTerm; 
                biasUpdate = sum(gc,2);
            end
        else
            throw(CORAerror('CORA:wrongFieldValue','options.nn.train.zonotope_weight_update',...
               "Only supported values are 'center' and 'extreme'!"));
        end
        % update weights and bias
        obj.updateGrad('W', weightsUpdate, options);
        obj.updateGrad('b', biasUpdate, options);

        % linear map of the out-going gradient
        if ~options.nn.interval_center
    	    gc = obj.W'*gc;
        end
        gG = pagemtimes(obj.W',gG);

        % Clear backprop storage.
        clear 'obj.backprop.store';
    end
end

% Auxiliary functions -----------------------------------------------------

methods
    function names = getLearnableParamNames(obj)
        % list of learnable properties
        names = {'W', 'b'};
    end

    function [mask,keepIdx,dropIdx] = getDropMask(obj,x,dropFactor)
        % Get the size of the input.
        [n,batchSize] = size(x);
        % Compute random permutation of the dimensions for each element in 
        % the batch.
        [~,dDims] = sort(rand([n batchSize],'like',x),1);
        % Convert to linear indices.
        dDimsIdx = reshape(sub2ind(size(x),dDims, ...
            repmat(1:batchSize,n,1)),size(dDims));
        % Get drop factor.
        % dropFactor = obj.dropFactor;
        % Compute number of dimensions to keep.
        numDimsKeep = ceil(n*dropFactor);
        % Obtain dimensions to keep and which to set to 0.
        keepIdx = dDimsIdx(1:numDimsKeep,:);
        dropIdx = dDimsIdx(numDimsKeep+1:end,:);
        % Construct the mask.
        mask = zeros([n batchSize],'like',x);
        % Scale remaining dimensions s.t. their sum remains constant.
        mask(keepIdx) = 1/(1 - dropFactor);
    end
end

% Auxiliary functions -----------------------------------------------------

methods

    function outputSize = getOutputSize(obj, inputSize)
        outputSize = [size(obj.W, 1), 1];
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = size(obj.W, 2);
        nout = size(obj.W, 1);
    end

    function fieldStruct = getFieldStruct(obj)
        fieldStruct = struct;
        fieldStruct.size_W = size(obj.W); % is lost for vectors in json
        fieldStruct.W = obj.W;
        fieldStruct.b = obj.b;
        fieldStruct.d = obj.d;
    end

end

end

% ------------------------------ END OF CODE ------------------------------
