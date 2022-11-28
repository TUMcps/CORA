classdef (Abstract) nnLayer < handle
% nnLayer - abstract class for nn layers
%
% Syntax:
%    abstract class
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

properties (Abstract, Constant)
    type % implemented in subclasses
end
properties
    name % could be helpful for debugging
end

methods
    function obj = nnLayer(name)
        if nargin < 1
            name = nnLayer.type;
        end

        obj.name = name;
    end
end

methods (Abstract)
    [nin, nout] = getNumNeurons(obj)
    % evaluation functions for each set representation
    evaluateNumeric(obj, input)
    evaluateZonotope(obj, Z, evParams)
    evaluatePolyZonotope(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
    evaluateTaylm(obj, input, evParams)
    evaluateConZonotope(obj, c, G, C, d, l, u, options, evParams)
end
end

%------------- END OF CODE --------------