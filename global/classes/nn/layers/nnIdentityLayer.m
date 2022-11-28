classdef nnIdentityLayer < nnLayer
% nnIdentityLayer - class for identity layer
%    This layer is usually not necessary but can be helpful for e.g.
%    conversion from neuralNetworkOld model or keeping order in layers
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

% Author:       Tobias Ladner
% Written:      28-March-2022
% Last update:  ---
% Last revision:10-August-2022 (renamed)

%------------- BEGIN CODE --------------

properties (Constant)
    type = "IdentityLayer"
end

methods
    % constructor
    function obj = nnIdentityLayer(name)
        if nargin < 1
            name = nnIdentityLayer.type;
        end
        % call super class constructor
        obj@nnLayer(name)
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = [];
        nout = [];
    end

    % evaluate
    function input = evaluateNumeric(obj, input)
        % return identity
    end

    function Z = evaluateZonotope(obj, Z, evParams)
        % return identity
    end

    function [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
        % return identity
    end

    function input = evaluateTaylm(obj, input, evParams)
        % return identity
    end

    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, j, options, evParams)
        % return identity
    end
end
end

%------------- END OF CODE --------------