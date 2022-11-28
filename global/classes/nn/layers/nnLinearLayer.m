classdef nnLinearLayer < nnLayer
% nnLinearLayer - class for linear layers
%
% Syntax:
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

% Author:       Tobias Ladner
% Written:      28-March-2022
%               23-November-2022 (polish)
% Last update:  ---
% Last revision:10-August-2022 (renamed)

%------------- BEGIN CODE --------------

properties (Constant)
    type = "LinearLayer"
end

properties
    W, b
end

methods
    % constructor
    function obj = nnLinearLayer(W, b, name)
        if nargin < 2
            throw(CORAerror('CORA:notEnoughInputArgs', 2));
        end

        if ~(isa(W, 'double') && isa(b, 'double'))
           throw(CORAerror('CORA:wrongInputInConstructor', ...
               'W and b should be of type double.'));
        end
        if ~all(size(b) == [size(W, 1), 1])
           throw(CORAerror('CORA:wrongInputInConstructor', ...
               'The dimensions of W and b should match.'));
        end

        if nargin < 3
            name = nnLinearLayer.type;
        end
        % call super class constructor
        obj@nnLayer(name)

        obj.W = W;
        obj.b = b;
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = size(obj.W, 2);
        nout = size(obj.W, 1);
    end

    % evaluate
    function r = evaluateNumeric(obj, input)
        r = obj.W * input + obj.b;
    end

    function r = evaluateZonotope(obj, Z, evParams)
        Z = obj.W * Z;
        Z(:, 1) = Z(:, 1) + obj.b;
        r = Z;
    end

    function [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
        c = obj.W * c + obj.b;
        G = obj.W * G;
        Grest = obj.W * Grest;
    end

    function r = evaluateTaylm(obj, input, evParams)
        r = obj.W * input + obj.b;
    end

    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options, evParams)
        c = obj.W * c + obj.b;
        G = obj.W * G;
    end
end
end

%------------- END OF CODE --------------