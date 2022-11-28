classdef nnSigmoidLayer < nnSShapeLayer
% nnSigmoidLayer - class for Sigmoid layers
%
% Syntax:
%    obj = nnSigmoidLayer(name)
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
    type = "SigmoidLayer"
end

methods
    % constructor
    function obj = nnSigmoidLayer(name)
        if nargin < 1
            name = nnSigmoidLayer.type;
        end
        % call super class constructor
        obj@nnSShapeLayer(name)
    end

    % evaluate
    function r = evaluateNumeric(obj, input)
        % use tanh for numeric stability
        r = tanh(input/2) / 2 + 0.5;
    end
end
end

%------------- END OF CODE --------------