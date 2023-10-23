classdef nnExpLayer < nnActivationLayer
% nnExpLayer - class for exp layers
%
% Syntax:
%    obj = nnExpLayer(name)
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

% Authors:       Tobias Ladner
% Written:       19-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

methods
    % constructor
    function obj = nnExpLayer(name)
        if nargin < 2
            name = [];
        end
        % call super class constructor
        obj@nnActivationLayer(name)
    end

    function der1 = getDerBounds(obj, l, u)
        % df_l and df_u as lower and upper bound for the derivative
        der1 = interval(exp(l), exp(u));
    end
end

% evaluate ----------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})
    % numeric
    function [r, obj] = evaluateNumeric(obj, input, evParams)
        r = exp(input);
    end
end

methods (Access=protected)
    function [coeffs, d] = computeApproxPolyCustom(obj, l, u, order, poly_method)
        % implement custom polynomial computation in subclass
        coeffs = []; d = [];
        
        f = obj.f;
        df = obj.getDf(1);

        if strcmp(poly_method, 'singh')
            if order == 1
                % according to [1, Theorem 3.2] where we interpret
                % exp as lower part of an S-curve
                lambda = min(df(l), df(u));
                mu1 = 0.5 * (f(u) + f(l) - lambda * (u + l));
                mu2 = 0.5 * (f(u) - f(l) - lambda * (u - l));
                coeffs = [lambda, mu1];
                d = mu2;
            end
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
