classdef nnInvSqrtRootLayer < nnActivationLayer
% nnInvSqrtRootLayer - class for inverse square root layers
%    only defined on x\in(0,Inf]
%
% Syntax:
%    obj = nnInvSqrtRootLayer(name)
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
% Written:       15-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

methods
    % constructor
    function obj = nnInvSqrtRootLayer(name)
        if nargin < 1
            name = [];
        end
        % call super class constructor
        obj@nnActivationLayer(name)
    end

    function df_i = getDf(obj, i)
        if i == 0
            df_i = obj.f;
        elseif i == 1
            df_i = @(x) -1 / 2 .* x.^(-3/2);
        else
            throw(CORAerror('CORA:notSupported','Higher order derivative than first derivative.'))
        end
    end

    function [df_l, df_u] = getDerBounds(obj, l, u)
        df_l = obj.df(l);
        df_u = obj.df(u);
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})
    % numeric
    function r = evaluateNumeric(obj, input, options)
        r = input.^(-1/2);
    end
end

methods (Access=protected)
    function [coeffs, d] = computeApproxPolyCustom(obj, l, u, order, poly_method)
        % implement custom polynomial computation in subclass
        
        % check domain
        if l <= 0
            throw(CORAerror('CORA:outOfDomain','validDomain','x\\in(0,Inf]'))
        end

        % init
        coeffs = []; d = [];
        f = obj.f;
        df = obj.getDf(1);

        if strcmp(poly_method, 'singh')
            if order == 1
                % analogous to [1, Theorem 3.2]
                lambda = df(u);
                mu1 = 0.5 * (f(u) + f(l) - lambda * (u + l));
                mu2 = 0.5 * (f(l) - (f(u) - lambda * (u - l)));
                coeffs = [lambda, mu1];
                d = mu2;
            end
        end
    end
end

methods
    function [coeffs, d] = computeApproxError(obj, l, u, coeffs)
        % compute order
        order = numel(coeffs)-1;

        % check domain
        if l <= 0
            throw(CORAerror('CORA:outOfDomain','validDomain','x\\in(0,Inf]'))
        end
        
        if order == 1
            % use analytic solution (either at {l,u,x0})
            x0 = (-1/(2*coeffs(1))).^(2/3);
            
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
            [coeffs, d] = computeApproxError@nnActivationLayer(obj,l,u,coeffs);
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
