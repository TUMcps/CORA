classdef nnElementwiseAffineLayer < nnLayer
% nnElementwiseAffineLayer - class for elementwise affine layers
%
% Syntax:
%    obj = nnElementwiseAffineLayer(scale)
%    obj = nnElementwiseAffineLayer(scale, offset)
%    obj = nnElementwiseAffineLayer(scale, offset, name)
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

% Authors:       Tobias Ladner, Lukas Koller
% Written:       30-March-2022
% Last update:   14-December-2022 (variable input tests, inputArgsCheck)
%                21-March-2024 (batchZonotope for training)
% Last revision: 10-August-2022 (renamed)

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    scale, offset
end

methods
    % constructor
    function obj = nnElementwiseAffineLayer(varargin)
        % parse input
        [scale, offset, name] = setDefaultValues({1, 0, []}, varargin);
        inputArgsCheck({ ...
            {scale, 'att', 'numeric'}
            {offset, 'att', 'numeric'}
        });

        % check dims
        if size(scale, 2) > 1 || size(offset, 2) > 1
           throw(CORAerror('CORA:wrongInputInConstructor', ...
               'Scale and offset should be column vectors.'));
        end
        if length(scale) > 1 && length(offset) > 1 && ...
            length(scale) ~= length(offset)
           throw(CORAerror('CORA:wrongInputInConstructor', ...
               'The dimensions of scale and offset should match or be scalar values.'));
        end

        % call super class constructor
        obj@nnLayer(name)

        obj.scale = double(scale);
        obj.offset = double(offset);
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = [];
        nout = [];
    end

    function outputSize = getOutputSize(obj, inputSize)
        outputSize = inputSize;
    end

    function castWeights(obj, x)
        % Callback when data type of learnable parameters was changed
        obj.scale = cast(obj.scale,'like',x);
        obj.offset = cast(obj.offset,'like',x);
    end
end

% evaluate ----------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})

    % numeric
    function r = evaluateNumeric(obj, input, options)
        [scale,offset] = obj.aux_getScaleAndOffset();
        r = scale(:) .* input + offset(:);
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, options)
        [scale,offset] = obj.aux_getScaleAndOffset();
        S = scale(:)' .* S;

        if options.nn.store_sensitivity
            % Store the gradient (used for the sensitivity computation).
            obj.sensitivity = S;
        end
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
        [scale,offset] = obj.aux_getScaleAndOffset();
        c = scale(:) .* c + offset;
        G = scale(:) .* G;
        GI = scale(:) .* GI;
    end

    % interval 
    function bounds = evaluateInterval(obj, bounds, options)
        [scale,offset] = obj.aux_getScaleAndOffset();
        l_ = scale.*bounds.inf + offset;
        u_ = scale.*bounds.sup + offset;
        bounds = interval(min(l_,u_),max(l_,u_));
    end

    % zonotope batch (for training)
    function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
        [scale,offset] = obj.aux_getScaleAndOffset();
        % Add the offset.
        c = scale(:).*c + offset(:);
        if options.nn.interval_center
            % Flip bounds in case the scale is negative.
            c = [min(c,[],2) max(c,[],2)];
        end
        % Scale the generators.
        G = scale(:).*G;
    end

    % taylm
    function r = evaluateTaylm(obj, input, options)
        [scale,offset] = obj.aux_getScaleAndOffset();
        r = scale * input + offset;
    end

    % conZonotope
    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options)
        [scale,offset] = obj.aux_getScaleAndOffset();
        c = scale * c + offset;
        G = scale * G;
    end

    % backprop ------------------------------------------------------------

    % numeric
    function grad_in = backpropNumeric(obj, input, grad_out, options, updateWeights)
        [scale,offset] = obj.aux_getScaleAndOffset();
        grad_in = scale .* grad_out;
    end

    % interval batch
    function [gl, gu] = backpropIntervalBatch(obj, l, u, gl, gu, options, updateWeights)
        [scale,offset] = obj.aux_getScaleAndOffset();
        gl = scale.*gl;
        gu = scale.*gu;
    end
    
    % zonotope batch
    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options, updateWeights)
        [scale,offset] = obj.aux_getScaleAndOffset();
        gc = scale.*gc;
        gG = scale.*gG;
    end
end


% Auxiliary functions -----------------------------------------------------

methods
    function fieldStruct = getFieldStruct(obj)
        fieldStruct = struct;
        fieldStruct.scale = obj.scale;
        fieldStruct.offset = obj.offset;
    end
end

% protected methods
methods (Access = protected)

    % Pad scale or offset channel-wise.
    function p = aux_getPaddedParameter(obj, p, varargin)
        [inImgSize] = setDefaultValues({obj.inputSize}, varargin);
        % Ensure the image size has at least 3 dimensions.
        ndims = length(inImgSize);
        if ndims < 3
            inImgSize = [inImgSize ones(3 - ndims)];
        end

        % Compute number of spacial dimensions in the feature map.
        spacDim = prod(inImgSize(1:2));
        % Obtain number of output channels.
        out_c = inImgSize(3);

        if numel(p) < spacDim*out_c
            if isempty(p)
                p = zeros(out_c,1,'like',p);
            elseif isscalar(p)
                p = repmat(p(:),[out_c 1]);
            end
    
            % Expand the parameter vector to output size.
            p = repelem(p(:), spacDim, 1);
        else
            p = p(:);
        end
    end

    function [scale,offset] = aux_getScaleAndOffset(obj)
        if isfield(obj.backprop.store,'scale') ...
                && isfield(obj.backprop.store,'offset')
            % Return stored parameters.
            scale = obj.backprop.store.scale;
            offset = obj.backprop.store.offset;
        else
            % Obtain input size.
            inImgSize = obj.inputSize;
            % Pad scale and offset to match input size; i.e., channel-wise
            % scale and offset.
            scale = aux_getPaddedParameter(obj, obj.scale, inImgSize);
            offset = aux_getPaddedParameter(obj, obj.offset, inImgSize);
            % Store the parameters.
            obj.backprop.store.scale = scale;
            obj.backprop.store.offset = offset;
        end

    end
end

end

% ------------------------------ END OF CODE ------------------------------
