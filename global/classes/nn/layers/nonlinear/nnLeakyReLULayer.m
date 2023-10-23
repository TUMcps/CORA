classdef nnLeakyReLULayer < nnActivationLayer
% nnLeakyReLULayer - class for LeakyReLU layers
%
% Syntax:
%    obj = nnLeakyReLULayer(alpha, name)
%
% Inputs:
%    alpha - slope of the LeakyReLU for x<0, defaults to 0.01
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Singh, G., et al. "Fast and Effective Robustness Certification"
%    [2] Kochdumper, Niklas, et al. "Open-and closed-loop neural network 
%       verification using polynomial zonotopes." 
%       arXiv preprint arXiv:2207.02715 (2022).
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Sebastian Sigl, Tobias Ladner
% Written:       11-June-2022
% Last update:   16-February-2023 (TL, combined approx_type)
% Last revision: 10-August-2022 (renamed)
%                26-May-2023

% ------------------------------ BEGIN CODE -------------------------------

properties
    alpha
end

methods
    % constructor
    function obj = nnLeakyReLULayer(alpha, name)
        if nargin < 1
            alpha = 0.01;
        end
        if nargin < 2
            name = [];
        end
        
        inputArgsCheck({{alpha,'att','numeric','scalar'}});

        % call super class constructor
        obj@nnActivationLayer(name)
        obj.alpha = alpha;
    end
end

% evaluate ----------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})
    % numeric
    function r = evaluateNumeric(obj, input, evParams)
        r = max(obj.alpha.*input, input);
    end
end

% Auxiliary functions -----------------------------------------------------

methods
    function df_i = getDf(obj, i)
        if i == 0
            df_i = obj.f;
        elseif i == 1
            df_i = @(x) 1 * (x > 0) + obj.alpha * (x <= 0);
        else
            df_i = @(x) 0;
        end
    end

    function [df_l,df_u] = getDerBounds(obj, l, u)
        % df_l and df_u as lower and upper bound for the derivative
        % case distinction for l
        if l <= 0
            df_l = obj.alpha;
        else
            df_l = 1;
        end

        % case distinction for u
        if u < 0
            df_u = obj.alpha;
        else
            df_u = 1;
        end
    end

    function [coeffs, d] = computeApproxPoly(obj, l, u, varargin)
        % computes an approximating polynomial and respective error bound
        % exploit piecewise linearity of nnLeakyReLULayer
  
        % check if ReLU can be computed exactly
        if u <= 0
            coeffs = [obj.alpha, 0];
            d = 0; % no approximation error!

        elseif 0 <= l
            % identity
            coeffs = [1, 0];
            d = 0; % no approximation error!

        else % l < 0 < u
            [coeffs, d] = computeApproxPoly@nnActivationLayer(obj, l, u, varargin{:});
        end
    end

    function [coeffs, d] = computeApproxError(obj, l, u, coeffs)
        % error can be computed exactly by checking each linear part
        % according to [2, Sec. 3.2]

        % x < 0: p(x) - alpha*x
        [diffl1,diffu1] = minMaxDiffPoly(coeffs,[obj.alpha,0],l,0);
        
        % x > 0: p(x) - 1*x
        [diffl2,diffu2] = minMaxDiffPoly(coeffs,[1,0],0,u);
        
        % compute final approx error
        diffl = min(diffl1,diffl2);
        diffu = max(diffu1,diffu2);
        diffc = (diffu+diffl)/2;
        coeffs(end) = coeffs(end) - diffc;
        d = diffu-diffc; % error is radius then.
    end
end

methods(Access=protected)

    function [coeffs, d] = computeApproxPolyCustom(obj, l, u, order, poly_method)
        % implement custom polynomial computation in subclass
        coeffs = []; d = [];

        if strcmp(poly_method, 'singh')
            if order == 1
                % according to [1, Theorem 3.1]
                lambda = (u - obj.alpha * l) / (u - l);
                mu = 0.5 * (u - ((u - obj.alpha * l) / (u - l)) * u);
                coeffs = [lambda, mu];
                d = mu;
                return
            elseif order == 2
                % according to [2, Sec. 3.1]
                c_a = u * (1 - obj.alpha) / (u - l)^2;
                c_b = obj.alpha - 2 * c_a * l;
                c_c = c_a * l^2;
                coeffs = [c_a, c_b, c_c];
            end
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
