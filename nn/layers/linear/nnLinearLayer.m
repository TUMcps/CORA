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
        [b, name, areParamsLearnable] = ...
            setDefaultValues({0, [], true}, varargin);
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
        obj@nnLayer(name,areParamsLearnable)

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
    function S = evaluateSensitivity(obj, S, options)
        % Use pagemtimes to compute sensitivity simultaneously for an
        % entire batch.
        S = pagemtimes(S,obj.W);

        if options.nn.store_sensitivity
            % Store the gradient (used for the sensitivity computation).
            obj.sensitivity = S;
        end
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
            % The center consits of a lower and an upper bound; propagte
            % the center as an interval.
            cl = reshape(c(:,1,:),[n batchSize]);
            cu = reshape(c(:,2,:),[n batchSize]);
            c = obj.evaluateInterval(interval(min(cl,cu),max(cl,cu)));
            % Combine the center bounds.
            c = permute(cat(3,c.inf,c.sup),[1 3 2]);
        else
            c = obj.W*c + obj.b;
        end
        % Propagate the generator matrix.
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

    function grad_in = backpropNumeric(obj, input, grad_out, options, updateWeights)
        if updateWeights
            % update weights and bias
            obj.updateGrad('W', grad_out * input', options);
            obj.updateGrad('b', sum(grad_out, 2), options);
        end
        % backprop gradient
        grad_in = obj.W' * grad_out;
    end

    function [gl, gu] = backpropIntervalBatch(obj, l, u, gl, gu, options, updateWeights)
        % see Gowal et al. 2019
        mu = (u + l)/2;
        r = (u - l)/2;

        if updateWeights
            % update weights and bias
            obj.updateGrad('W', (gu + gl)*mu' ...
                + (gu - gl)*r'.*sign(obj.W), options);
            obj.updateGrad('b', sum(gu + gl, 2), options);
        end

        % backprop gradient
        dmu = obj.W' * (gu + gl)/2;
        dr = abs(obj.W') * (gu - gl)/2;
        gl = dmu - dr;
        gu = dmu + dr;
    end

    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options, updateWeights)
        [n,numGen,batchSize] = size(G);
        % obtain indices of active generators
        genIds = obj.genIds;

        if options.nn.interval_center
            % Compute gradient of the interval center.
            cl = reshape(c(:,1,:),[n batchSize]);
            cu = reshape(c(:,2,:),[n batchSize]);
            gl = reshape(gc(:,1,:),[size(gc,1) batchSize]);
            gu = reshape(gc(:,2,:),[size(gc,1) batchSize]);
            % Compute the out going gradient of the center.
            [gl, gu] = backpropIntervalBatch(obj, cl, cu, gl, gu, ...
                options, updateWeights);
            % Construct the out going interval center gradient.
            gc = permute(cat(3,gl,gu),[1 3 2]);
            % The center term and bias updated are already applied.
            centerTerm = 0;
            biasUpdate = 0;
        else
            % Compute outer product between centers.
            centerTerm = gc*c';
            % Compute bias update.
            biasUpdate = sum(gc,2);
            % Compute the outgoing gradient of the center.
            gc = obj.W'*gc;
        end

        if strcmp(options.nn.train.zonotope_weight_update,'center')
            % Use the center to update the weights and biases. There is no
            % generator term.
            gensTerm = 0;
        elseif strcmp(options.nn.train.zonotope_weight_update,'sum')
            % Compute the outer product between generator matrices.
            gensTerm = sum(pagemtimes(gG(:,genIds,:),'none', ...
                G(:,genIds,:),'transpose'),3);
        else
            throw(CORAerror('CORA:wrongFieldValue', ...
                'options.nn.train.zonotope_weight_update',...
               "Only supported values are 'center' and 'sum'!"));
        end

        if updateWeights
            % Compute weights and bias update.
            weightsUpdate = centerTerm + gensTerm; 
            % Update weights and bias.
            obj.updateGrad('W', weightsUpdate, options);
            obj.updateGrad('b', biasUpdate, options);
        end

        % Apply the linear map of the out-going gradient of the generator
        % matrix.
        gG = pagemtimes(obj.W',gG);
    end
end

% Auxiliary functions -----------------------------------------------------

methods
    function names = getParamNames(obj)
        % List parameters.
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
