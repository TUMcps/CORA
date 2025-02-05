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
    function S = evaluateSensitivity(obj, S, x, options)
        [scale,offset] = obj.aux_getScaleAndOffset();
        S = scale(:)' .* S;
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, optionss)
        c = obj.scale * c + obj.offset;
        G = obj.scale * G;
        GI = obj.scale * GI;
    end

    % interval 
    function bounds = evaluateInterval(obj, bounds, options)
        l = obj.scale.*bounds.inf + obj.offset;
        u = obj.scale.*bounds.sup + obj.offset;
        bounds = interval(l,u);
    end

    % zonotope batch (for training)
    function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
        [scale,offset] = obj.aux_getScaleAndOffset();
        if options.nn.interval_center
            % Flip bounds in case the scale is negative.
            mask = [(scale(:) < 0) ~(scale(:) < 0)];
            c_ = permute(c,[3 1 2]);
            c = permute(cat(3,c_(:,mask),c_(:,~mask)),[2 3 1]);
        end
        % Add the offset.
        c = scale(:).*c + offset(:);
        if options.nn.interval_center
            % Flip bounds in case the scale is negative.
            c = sort(c,2);
        end
        % Scale the generators.
        G = scale(:).*G;
    end

    % taylm
    function r = evaluateTaylm(obj, input, options)
        r = obj.scale * input + obj.offset;
    end

    % conZonotope
    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options)
        c = obj.scale * c + obj.offset;
        G = obj.scale * G;
    end

    % backprop ------------------------------------------------------------

    % numeric
    function grad_in = backpropNumeric(obj, input, grad_out, options)
        grad_in = obj.scale .* grad_out;
    end

    % interval batch
    function [gl, gu] = backpropIntervalBatch(obj, l, u, gl, gu, options)
        gl = obj.scale.*gl;
        gu = obj.scale.*gu;
    end

    % zonotope batch
    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options)
        gc = obj.scale.*gc;
        gG = obj.scale.*gG;
    end
end

% protected methods
methods (Access = protected)
    function [scale,offset] = aux_getScaleAndOffset(obj)
        % read out scale and offsets as vectors
        scale = obj.scale(:);
        offset = obj.offset(:);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
