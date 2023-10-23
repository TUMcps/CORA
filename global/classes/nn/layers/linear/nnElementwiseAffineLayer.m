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

% Authors:       Tobias Ladner
% Written:       30-March-2022
% Last update:   14-December-2022 (variable input tests, inputArgsCheck)
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

end

% evaluate ----------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})

    % numeric
    function r = evaluateNumeric(obj, input, evParams)
        r = obj.scale .* input + obj.offset;
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, evParams)
        S = obj.scale .* S;
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
        c = obj.scale * c + obj.offset;
        G = obj.scale * G;
        GI = obj.scale * GI;
    end

    % taylm
    function r = evaluateTaylm(obj, input, evParams)
        r = obj.scale * input + obj.offset;
    end

    % conZonotope
    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options, evParams)
        c = obj.scale * c + obj.offset;
        G = obj.scale * G;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
