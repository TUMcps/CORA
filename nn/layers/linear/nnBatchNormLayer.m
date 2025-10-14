classdef nnBatchNormLayer < nnElementwiseAffineLayer
% nnBatchNormLayer - class for batch normalization layers
%
% Syntax:
%    obj = nnBatchNormLayer(scale)
%    obj = nnBatchNormLayer(scale, offset)
%    obj = nnBatchNormLayer(scale, offset, name)
%
% Inputs:
%    scale - elementwise scale (scalar or matching dimension)
%    offset - elementwise offset (scalar or matching dimension)
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

% Authors:       Lukas Koller
% Written:       26-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    epsilon % small constant for numerical stability
    momentum % momentum for moving mean and variance
    movMean % trained mean
    movVar % trained variance
end

methods
    % constructor
    function obj = nnBatchNormLayer(varargin)
        % parse input
        [scale, offset, movVar, movMean, name, epsilon, momentum] = ...
            setDefaultValues({1, 0, 1, 0, [], 0.001, 0.99}, varargin);

        % call super class constructor
        obj@nnElementwiseAffineLayer(scale, offset, name)

        obj.epsilon = epsilon;
        obj.momentum = momentum;
        obj.movVar = movVar;
        obj.movMean = movMean;
    end

    function castWeights(obj, x)
        obj.castWeights@nnElementwiseAffineLayer(x)
        % Callback when data type of learnable parameters was changed
        obj.movMean = cast(obj.movMean,'like',x);
        obj.movVar = cast(obj.movVar,'like',x);
    end

    function names = getParamNames(obj)
        % list of learnable properties
        names = {'scale', 'offset'};
    end
end

% evaluate ----------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})

    % numeric
    function r = evaluateNumeric(obj, input, options)
        if (options.nn.train.backprop && ~options.nn.batch_norm_moving_stats) ...
                || options.nn.batch_norm_calc_stats
            % Obtain batch size.
            bs = size(input,2);
            % Obtain input size and normalization dimensions.
            [reshapeSize,normDims] = obj.aux_getNormDims(bs);
            % Compute mean and standard deviation across the batch.
            batchMean = mean(reshape(input,reshapeSize),normDims);
            batchVar = mean((reshape(input,reshapeSize) - batchMean).^2,normDims);
            % Update statistics.
            obj.movVar = obj.movVar * obj.momentum ...
                + batchVar * (1 - obj.momentum);
            obj.movMean = obj.movMean * obj.momentum ...
                + batchMean * (1 - obj.momentum);
        else
            % We are in inference-mode: Normalize the input with trained 
            % mean and variance.
            batchMean = obj.movMean;
            batchVar = obj.movVar;
        end
        % Correct pad and reshape the computed mean and variance.
        batchMean = obj.aux_getPaddedParameter(batchMean(:));
        batchVar = obj.aux_getPaddedParameter(batchVar(:));
        % Normalize the input with current mean and variance.
        isqrtVar = 1./sqrt(batchVar + obj.epsilon);
        input = (input - batchMean).*isqrtVar;
        % Apply scale and offset.
        r = obj.evaluateNumeric@nnElementwiseAffineLayer(input, options);

        if options.nn.train.backprop
            % Store the batch-normed input.
            obj.backprop.store.input_normed = input;
            obj.backprop.store.isqrtVar = isqrtVar;
        end
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, options)
        % Obtain the stored input.
        x = obj.backprop.store.input;

        % Apply learned scaling; use previously stored statistics.
        batchVar = obj.movVar; % obj.backprop.store.batchVar;
        batchVar = obj.aux_getPaddedParameter(batchVar);
        S = 1./sqrt(batchVar(:) + obj.epsilon)' .* S;
        % Apply scale and offset.
        S = obj.evaluateSensitivity@nnElementwiseAffineLayer(S, options);

        if options.nn.store_sensitivity
            % Store the gradient (used for the sensitivity computation).
            obj.sensitivity = S;
        end
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
        % Apply scale and offset.
        [c, G, GI, E, id, id_, ind, ind_] = ...
            obj.evaluatePolyZonotope@nnElementwiseAffineLayer( ...
                c, G, GI, E, id, id_, ind, ind_, options);
    end

    % interval 
    function bounds = evaluateInterval(obj, bounds, options)
        if (options.nn.train.backprop && ~options.nn.batch_norm_moving_stats) ...
                || options.nn.batch_norm_calc_stats
            % % We use the previously stored mean and variance of the nominal input.
            % batchMean = obj.backprop.store.batchMean;
            % batchVar = obj.backprop.store.batchVar;

            % Use the statistics of the center.
            c = 1/2*(bounds.sup + bounds.inf);
            % Obtain batch size.
            bs = size(c,2);
            % Obtain input size and normalization dimensions.
            [reshapeSize,normDims] = obj.aux_getNormDims(bs);
            % Compute mean and standard deviation across the batch.
            batchMean = mean(reshape(c,reshapeSize),normDims);
            batchVar = mean((reshape(c,reshapeSize) - batchMean).^2,normDims);
            % Update statistics.
            obj.movVar = obj.movVar * obj.momentum ...
                + batchVar * (1 - obj.momentum);
            obj.movMean = obj.movMean * obj.momentum ...
                + batchMean * (1 - obj.momentum);
        else
            % We are in inference-mode: Normalize the input with trained 
            % mean and variance.
            batchMean = obj.movMean;
            batchVar = obj.movVar;
        end
        % Correct pad and reshape the computed mean and variance.
        batchMean = obj.aux_getPaddedParameter(batchMean(:));
        batchVar = obj.aux_getPaddedParameter(batchVar(:));
        % Normalize the input with current mean and variance.
        isqrtVar = 1./sqrt(batchVar + obj.epsilon);
        bounds = interval( ...
            (bounds.inf - batchMean).*isqrtVar,...
            (bounds.sup - batchMean).*isqrtVar ...
        );

        % Apply scale and offset.
        bounds = obj.evaluateInterval@nnElementwiseAffineLayer(bounds, options);

        if options.nn.train.backprop
            % Store the batch-normed input.
            obj.backprop.store.bounds_normed = bounds;
            obj.backprop.store.isqrtVar = isqrtVar;
        end
    end

    % zonotope batch (for training)
    function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
        if (options.nn.train.backprop && ~options.nn.batch_norm_moving_stats) ...
                || options.nn.batch_norm_calc_stats
            if options.nn.interval_center
                c_ = reshape(1/2*sum(c,2),size(c,[1 3]));
            else
                c_ = c;
            end
            % Obtain batch size.
            bs = size(c_,2);
            % Obtain input size and normalization dimensions.
            [reshapeSize,normDims] = obj.aux_getNormDims(bs);
            % Compute mean and standard deviation across the batch.
            batchMean = mean(reshape(c_,reshapeSize),normDims);
            batchVar = mean((reshape(c_,reshapeSize) - batchMean).^2,normDims);
            % Update statistics.
            obj.movVar = obj.movVar * obj.momentum ...
                + batchVar * (1 - obj.momentum);
            obj.movMean = obj.movMean * obj.momentum ...
                + batchMean * (1 - obj.momentum);
        else
            % We are in inference-mode: Normalize the input with trained 
            % mean and variance.
            batchMean = obj.movMean;
            batchVar = obj.movVar;
        end
        % Correct pad and reshape the computed mean and variance.
        batchMean = obj.aux_getPaddedParameter(batchMean(:));
        batchVar = obj.aux_getPaddedParameter(batchVar(:));
        % Normalize the input with current mean and variance.
        isqrtVar = 1./sqrt(batchVar + obj.epsilon);
        c = (c - batchMean).*isqrtVar;
        G = G.*isqrtVar;

        % Apply scale and offset.
        [c, G] = obj.evaluateZonotopeBatch@nnElementwiseAffineLayer(c, G, options);

        if options.nn.train.backprop
            % Store the batch-normed input.
            obj.backprop.store.c_normed = c;
            obj.backprop.store.G_normed = G;
            obj.backprop.store.isqrtVar = isqrtVar;
        end
    end

    % backprop ------------------------------------------------------------

    function grad_in = backpropNumeric(obj, input, grad_out, options, updateWeights)
        % Obtain batch size.
        [~,bs] = size(input);

        % Obtain the stored batch-normed input.
        input_normed = obj.backprop.store.input_normed;
        % Compute Hadamard product between gradient and input.
        gradInput_ = grad_out.*input_normed;

        if updateWeights
            % Obtain input size and normalization dimensions.
            [reshapeSize,normDims] = obj.aux_getNormDims(bs);
            % Update offset parameter.
            gradOffset = sum(reshape(grad_out,reshapeSize),normDims);
            obj.updateGrad('offset',gradOffset,options);
            % Update scale parameter.
            gradScale = sum(reshape(gradInput_,reshapeSize),normDims);
            obj.updateGrad('scale',gradScale,options);
        end

        % Backprop gradient through scale and offset.
        grad_in = obj.backpropNumeric@nnElementwiseAffineLayer( ...
            input, grad_out, options, updateWeights);
        % Obtain stored batch statistics.
        isqrtVar = obj.backprop.store.isqrtVar;
        % Backprop through batch normalization.
        % grad_in = isqrtVar.*(grad_in - 1/bs*sum(grad_in,2) ...
        %     -1/bs*sum(gradInput_,2).*(input_normed - 1/bs*sum(input_normed,2)));
        grad_in = isqrtVar.*(grad_in - 1/bs*sum(grad_in,2) ...
            + 1/(2*bs^2)*isqrtVar.*sum(gradInput_,2));
    end

    function [gl, gu] = backpropIntervalBatch(obj, l, u, gl, gu, options, updateWeights)
        % Obtain batch size.
        [~,bs] = size(l);

        % Obtain the stored batch-normed input.
        bounds_normed = obj.backprop.store.bounds_normed;
        % Compute Hadamard product between gradient and input.
        glInput_ = gl.*bounds_normed.inf;
        guInput_ = gu.*bounds_normed.sup;

        if updateWeights
            % Obtain input size and normalization dimensions.
            [reshapeSize,normDims] = obj.aux_getNormDims(bs);
            % Update offset parameter.
            gradOffset = sum(reshape(gu + gl,reshapeSize),normDims);
            obj.updateGrad('offset',gradOffset,options);
            % Update scale parameter.
            gradScale = sum(reshape(glInput_ + guInput_,reshapeSize),normDims);
            obj.updateGrad('scale',gradScale,options);
        end

        % Backprop gradient through scale and offset.
        [gl, gu] = obj.backpropIntervalBatch@nnElementwiseAffineLayer( ...
            l, u, gl, gu, options, updateWeights);
        % Obtain stored batch statistics.
        isqrtVar = obj.backprop.store.isqrtVar;        
        % % Backprop through batch normalization; the statistics are computed 
        % % based on the nominal input, thus derivate of the mean and
        % % variance are 0.
        % gl = isqrtVar.*gl;
        % gu = isqrtVar.*gu;

        % grad_c = 1/2*(gu + gl);
        % grad_c = isqrtVar.*(grad_c - 1/bs*sum(grad_c,2) ...
        %     -1/bs*sum(gradInput_,2).*(input_normed - 1/bs*sum(input_normed,2)));

        sumg = 1/2*(1/bs*sum(gl + gu,2) ...
            - 1/(2*bs^2)*isqrtVar.*sum(glInput_ + guInput_,2));
        gl = isqrtVar.*(gl - sumg);
        gu = isqrtVar.*(gu - sumg);
    end

    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options, updateWeights)
        % Obtain batch size.
        [~,~,bs] = size(G);
        % Obtain the indices of the relevant generators.
        genIds = obj.genIds;

        % Obtain the stored batch-normed input.
        c_normed = obj.backprop.store.c_normed;
        G_normed = obj.backprop.store.G_normed(:,genIds,:);
        % Compute Hadamard product between gradient and input.
        if options.nn.interval_center
            gradc_ = reshape(sum(gc.*c_normed,2),size(c,[1 3]));
        else
            gradc_ = gc.*c_normed;
        end
        gradG_ = reshape(sum(gG(:,genIds,:).*G_normed,2),size(c,[1 3]));

        if updateWeights
            % Obtain input size and normalization dimensions.
            [reshapeSize,normDims] = obj.aux_getNormDims(bs);
            % Update offset parameter.
            if options.nn.interval_center
                gc_ = reshape(sum(gc,2),size(c,[1 3]));
            else
                gc_ = gc;
            end
            gradOffset = sum(reshape(gc_,reshapeSize),normDims);
            obj.updateGrad('offset',gradOffset,options);
            % Update scale parameter.
            gradScale = sum(reshape(gradc_,reshapeSize),normDims);
            obj.updateGrad('scale',gradScale,options);
        end

        % Backprop gradient through scale and offset.
        [gc, gG] = obj.backpropZonotopeBatch@nnElementwiseAffineLayer( ...
            c, G, gc, gG, options, updateWeights);
        % Obtain stored batch statistics.
        isqrtVar = obj.backprop.store.isqrtVar;
        % Backprop through batch normalization.
        % gc = isqrtVar.*(gc - 1/bs*sum(gc,2) ...
        %     -1/bs*sum(gradc_,2).*(c_normed - 1/bs*sum(c_normed,2)));
        % gG = isqrtVar.*gG;

        if options.nn.interval_center
            sumg = 1/2*(1/bs*sum(gc,[2 3]) ...
                - 1/(2*bs^2)*isqrtVar.*sum(gradc_ + gradG_,2));
        else
            sumg = (1/bs*sum(gc,2) ...
                + 1/(2*bs^2)*isqrtVar.*sum(gradc_ + gradG_,2));
        end
        gc = isqrtVar.*(gc - sumg);
        gG = isqrtVar.*gG;
    end
end

methods (Access = protected)

    function [reshapeSize,normDims] = aux_getNormDims(obj,bs)
        % Obtain input size.
        inImgSize = obj.inputSize;
        % Construct reshape size.
        reshapeSize = [inImgSize bs];
        % Obtain normalization axis.
        if numel(inImgSize) > 2
            % Normalize across channels.
            axis = 3;
        else
            % For vector input we normalize across the input
            % dimensions.
            axis = 1;
        end
        normDims = setdiff(1:length(inImgSize) + 1,axis);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
