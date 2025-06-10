classdef nnRootLayer < nnSShapeLayer
% nnRootLayer - class for square root layers
%
% Syntax:
%    obj = nnRootLayer(name)
%
% Inputs:
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner
% Written:       04-July-2022
% Last update:   02-May-2025 (TL, added analytic solution)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    reg_polys = [];
end

methods
    % constructor
    function obj = nnRootLayer(name)
        if nargin < 1
            name = [];
        end
        % call super class constructor
        obj@nnSShapeLayer(name)
    end

    function df_i = getDf(obj, i)
        if i == 0
            df_i = obj.f;
        else
            df_i1 = obj.getDf(i-1);
            df_i = @(x) 1 / 2 .* df_i1(x);
        end
    end

    function [df_l, df_u] = getDerBounds(obj, l, u)
        df_l = obj.df(u);
        df_u = obj.df(l);
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})
    % numeric
    function r = evaluateNumeric(obj, input, options)
        r = sqrt(input);
    end
end

% Auxiliary functions -----------------------------------------------------

methods
    function [coeffs, d] = computeApproxError(obj, l, u, coeffs)
        if numel(coeffs) == 2
            % use analytic solution (either at {l,u,x0})
            m = coeffs(1);
            x0 = 1/(2*m)^2;
            
            % sample possible solutions
            if l < x0 && x0 < u
                xs = [l,x0,u];
            else
                xs = [l,u];
            end

            % compute max error
            d_l = min(obj.f(xs) - (coeffs(1) .* xs + coeffs(2)));
            d_u = max(obj.f(xs) - (coeffs(1) .* xs + coeffs(2)));
            
            % do middle shift
            coeffs(end) = coeffs(end) + (d_l+d_u)/2;
            d = (d_u-d_l)/2;

        else
            [coeffs,d] = computeApproxError@nnSShapeLayer(obj,l,u,coeffs);
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
