classdef (Abstract) nnSShapeLayer < nnActivationLayer
% nnSShapeLayer - abstract class for S-shape like layers (sigmoid, tanh)
%
% Syntax:
%    nnSShapeLayer
%
% Inputs:
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Singh, G., et al. "Fast and Effective Robustness Certification"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   16-February-2023 (combined approx_type)
%                19-July-2023 (approx error: use sampling if very small)
% Last revision: 10-August-2022 (renamed)

% ------------------------------ BEGIN CODE -------------------------------

properties (Abstract)
    reg_polys % region polynomials
end

methods
    % constructor
    function obj = nnSShapeLayer(name)
        % call super class constructor
        obj@nnActivationLayer(name)
    end

    function [df_l,df_u] = getDerBounds(obj, l, u)
        df = obj.getDf(1);
        % df_l and df_u as lower and upper bound for the derivative
        % case distinction for all combinations of l and u
        if l <= 0 && u >= 0
            % global upper bound for the derivative of s-shaped functions
            df_u = df(0);
            % evaluate the lower bound in this case
            if abs(l) >= abs(u)
                df_l = df(l);
            else
                df_l = df(u);
            end
            % two leftover cases
        elseif l > 0 && u > 0
            df_l = df(u);
            df_u = df(l);
        elseif l < 0 && u < 0
            df_l = df(l);
            df_u = df(u);
        end
    
        % if l, u is very large, df is NaN
        if any(isnan([df_l, df_u]))
            df_l = 0;
            df_u = 0.0001;
        end
    end
end

% evaluate ----------------------------------------------------------------

% implemented in sup-/subclasses

% Auxiliary function ------------------------------------------------------

methods (Access=protected)
    function [coeffs, d] = computeApproxPolyCustom(obj, l, u, order, poly_method)
        % implement custom polynomial computation in subclass
        coeffs = []; d = [];
        
        f = obj.f;
        df = obj.getDf(1);

        if strcmp(poly_method, 'singh')
            if order == 1
                % according to [1, Theorem 3.2]
                lambda = min(df(l), df(u));
                mu1 = 0.5 * (f(u) + f(l) - lambda * (u + l));
                mu2 = 0.5 * (f(u) - f(l) - lambda * (u - l));
                coeffs = [lambda, mu1];
                d = mu2;
            end
        
        elseif strcmp(poly_method, "throw-catch")
            coeffs = nnHelper.calcAlternatingDerCoeffs(l, u, order, obj);
            
        end
    end
end

methods
    function [coeffs, d] = computeApproxError(obj, l, u, coeffs)
        % check if region polynomials are present
        P = obj.reg_polys;

        if isempty(P)
            % use sampling
            [coeffs,d] = computeApproxError@nnActivationLayer(obj,l,u,coeffs);

        else
            % use region polynomials
            P = P( ...
                ([P.l] <= l & l <= [P.u]) | ... % lower bound within region
                (l <= [P.l] & [P.u] <= u) | ... % region between bounds
                ([P.l] <= u & u <= [P.u])   ... % upper bound within region
            );
    
            % P should cover the entire domain
            if isempty(P)
                throw(CORAerror('CORA:emptySet',sprintf('Empty region polynomials for bounds=[%.4f,%.4f]. Make sure obj.polys covers the entire domain [-inf,inf].',l,u)))
            end
    
            % get approx error due to polys
            dP = max([P.d]);
    
            % find max difference between region polys and approximation poly
            diffl = Inf;
            diffu = -Inf;
            for i=1:length(P)
                [diffli,diffui] = minMaxDiffPoly( ...
                    coeffs, P(i).p, ...
                    max(l,P(i).l), ...
                    min(u,P(i).u));
                diffl = min(diffl,diffli);
                diffu = max(diffu,diffui);
            end
    
            % compute final approx error
            diffc = (diffl+diffu)/2;
            coeffs(end) = coeffs(end) - diffc;
            d = diffu-diffc + dP; % error is radius then.

            if diffu-diffc < 10*dP 
                % error is within order of magnitude of region poly. error
                % -> get tighter estimate using sampling
                [coeffs,d] = computeApproxError@nnActivationLayer(obj,l,u,coeffs);
            end
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
