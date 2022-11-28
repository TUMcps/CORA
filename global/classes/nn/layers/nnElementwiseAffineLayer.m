classdef nnElementwiseAffineLayer < nnLayer
% nnElementwiseAffineLayer - class for elementwise affine layers
%
% Syntax:
%    obj = nnElementwiseAffineLayer(scale, offset, name)
%
% Inputs:
%    scale - elementwise scale
%    offset - elementwise offset
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

% Author:       Tobias Ladner
% Written:      30-March-2022
% Last update:  ---
% Last revision:10-August-2022 (renamed)

%------------- BEGIN CODE --------------

properties (Constant)
    type = "ElementwiseAffineLayer"
end

properties
    scale, offset
end

methods
    % constructor
    function obj = nnElementwiseAffineLayer(scale, offset, name)
        if nargin < 3
            name = nnElementwiseAffineLayer.type;
        end
        % call super class constructor
        obj@nnLayer(name)

        obj.scale = scale;
        obj.offset = offset;
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = [];
        nout = [];
    end

    % evaluate
    function r = evaluateNumeric(obj, input)
        r = obj.scale .* input + obj.offset;
    end

    function r = evaluateZonotope(obj, Z, evParams)
        Z = obj.scale * Z;
        Z(:, 1) = Z(:, 1) + obj.offset;
        r = Z;
    end

    function [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
        c = obj.scale * c + obj.offset;
        G = obj.scale * G;
        Grest = obj.scale * Grest;
    end

    function r = evaluateTaylm(obj, input, evParams)
        r = obj.scale * input + obj.offset;
    end

    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options, evParams)
        c = obj.scale * c + obj.offset;
        G = obj.scale * G;
    end
end
end

%------------- END OF CODE --------------