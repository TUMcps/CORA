classdef nnTanhLayer < nnSShapeLayer
% nnTanhLayer - class for tanh layers
%
% Syntax:
%    obj = nnTanhLayer(name)
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
    type = "TanhLayer"
end

methods
    % constructor
    function obj = nnTanhLayer(name)
        if nargin < 1
            name = nnTanhLayer.type;
        end
        % call super class constructor
        obj@nnSShapeLayer(name)
    end

    function [r, obj] = evaluateNumeric(obj, input)
        r = tanh(input);
    end

end
end

%------------- END OF CODE --------------