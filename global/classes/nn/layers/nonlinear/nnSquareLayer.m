classdef nnSquareLayer < nnLayer
% nnSquareLayer - class for square layers
%
% Syntax:
%    obj = nnSquareLayer(name)
%
% Inputs:
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
% Written:       04-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    W, b
end

methods
    % constructor
    function obj = nnSquareLayer(name)
        if nargin < 1
            name = [];
        end
        % call super class constructor
        obj@nnLayer(name)
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = size(obj.W, 2);
        nout = size(obj.W, 1);
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})

    % numeric
    function r = evaluateNumeric(obj, input, evParams)
        r = input.^2;
    end

    % sensitvity
    function S = evaluateSensitivity(obj, S, x, evParams)
        S =  S * diag(2*x);
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
        [c, G, GI] = nnHelper.calcSquared(c, G, GI, E, c, G, GI, E, true);
        E = nnHelper.calcSquaredE(E, E, true);
        ind = find(prod(ones(size(E))-mod(E, 2), 1) == 1);
        ind_ = setdiff(1:size(E, 2), ind);
    end
end
end

% ------------------------------ END OF CODE ------------------------------
