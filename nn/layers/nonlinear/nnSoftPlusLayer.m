classdef nnSoftPlusLayer < nnActivationLayer
% nnSoftPlusLayer - class for SoftPlus layers, a smooth and continuous
%   approximation of ReLU.
%
% Syntax:
%    obj = nnSoftPlusLayer(name)
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
% See also: NeuralNetwork

% Authors:       Lukas Koller
% Written:       05-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    sharpness = 2       % sharpness parameter 
end

methods
    % constructor
    function obj = nnSoftPlusLayer(varargin)
        % Check number of input arguments.
        narginchk(0,2);
        % Parse input arguments.
        [sharpness,name] = setDefaultValues({1, []}, varargin);
        % call super class constructor
        obj@nnActivationLayer(name)
        obj.sharpness = sharpness;
    end

    function [df_l,df_u] = getDerBounds(obj, l, u)
        % Obtain the derivative.
        df = obj.getDf(1);
        % The derivative is monotonic.
        df_l = df(l);
        df_u = df(u);
    end
end

% evaluate ----------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})
    % Numeric evaluation.
    function [r, obj] = evaluateNumeric(obj, input, options)
        % Obtain the sharpness parameter.
        k = obj.sharpness;
        % Compute the output.
        r = (1/k).*log(1 + exp(k*input));
    end
end

methods (Access=protected)
    function [xs,dxsdm] = computeExtremePointsBatch(obj, m, options)
        % Obtain the sharpness parameter.
        k = obj.sharpness;
        % The slope has to be 0 < m < 1; specify an epsilon to avoid
        % numerical issues.
        tol = 1e-8;
        m = min(tol,max(m,1-tol));
        % Compute the exteme points.
        xs = (1/k).*log(-m./(m-1));
        % Compute the derivatie of the extreme points w.r.t. the slope.
        dxsdm = (1/k).*1./(m - m.^2);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
