classdef nnIdentityLayer < nnLayer
% nnIdentityLayer - class for identity layer
%    This layer is usually not necessary but can sometimes be helpful
%
% Syntax:
%    obj = nnIdentityLayer(name)
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
% Written:       28-March-2022
% Last update:   ---
% Last revision: 10-August-2022 (renamed)

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

methods
    % constructor
    function obj = nnIdentityLayer(name)
        if nargin < 1
            name = [];
        end
        % call super class constructor
        obj@nnLayer(name)
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
    function input = evaluateNumeric(obj, input, evParams)
        % return identity
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, evParams)
       % return identity
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
        % return identity
    end

    % taylm
    function input = evaluateTaylm(obj, input, evParams)
        % return identity
    end

    % conZonotope
    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, j, options, evParams)
        % return identity
    end
end

end

% ------------------------------ END OF CODE ------------------------------
