classdef (Abstract) nnSShapeLayer < nnActivationLayer
% nnSShapeLayer - abstract class for S-shape like layers (sigmoid, tanh)
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
% References:
%    [1] Singh, G., et al. "Fast and Effective Robustness Certification"
%    [2] R. Ivanov, et al. "Verisig 2.0: Verification of Neural Network
%        Controllers Using Taylor Model Preconditioning", 2021
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
end

properties
    sym_dfs = cell(0, 1)
end

methods
    % constructor
    function obj = nnSShapeLayer(name)
        % call super class constructor
        obj@nnActivationLayer(name)
    end

    function df_i = getDf(obj, i)
        if i == 0
            df_i = obj.f;
        elseif i <= size(obj.dfs, 1)
            df_i = obj.dfs{i};
        elseif i <= size(obj.sym_dfs, 1)
            df_i = matlabFunction(obj.sym_dfs{i});
        else
            syms x real
            if i == 1
                obj.sym_dfs{1, 1} = simplify(diff(obj.f(x)));
            end

            for j = 2:i
                obj.sym_dfs{j, 1} = simplify(diff(obj.sym_dfs{j-1, 1}));
            end

            df_i = matlabFunction(obj.sym_dfs{i});
        end
    end

    % evaluates
    function [res, d] = evaluateZonotopeNeuron(obj, input)
        % enclose the hyperbolic tangent activation function with a zonotope
        % according to Theorem 3.2 in [1]

        % hyperbolic tangent function and its derivative
        f = obj.f;
        df = obj.df;

        % compute lower and upper bound
        l = input(1) - sum(abs(input(2:end)));
        u = input(1) + sum(abs(input(2:end)));

        % compute auxiliary functions
        lambda = min(df(l), df(u));
        mu1 = 0.5 * (f(u) + f(l) - lambda * (u + l));
        mu2 = 0.5 * (f(u) - f(l) - lambda * (u - l));

        % compute output
        if l == u
            res = zeros(1, length(input));
            res(1) = f(u);
            d = 0;
        else
            res = lambda * input;
            res(1) = res(1) + mu1;
            d = mu2;
        end
    end

    function [c, G, Grest, d] = evaluatePolyZonotopeNeuronLin(obj, c, G, Grest, E, ind, ind_, approx)
        % enclose the hyperbolic tangent activation function with a polynomial
        % zonotope according to Theorem 3.2 in [1]

        % hyperbolic tangent function and its derivative
        f = obj.f;
        df = obj.df;

        % compute lower and upper bound
        [l, u] = nnHelper.compBoundsPolyZono(c, G, Grest, E, ind, ind_, approx);

        % compute auxiliary functions
        lambda = min(df(l), df(u));
        mu1 = 0.5 * (f(u) + f(l) - lambda * (u + l));
        mu2 = 0.5 * (f(u) - f(l) - lambda * (u - l));

        % compute output
        if l == u
            G = zeros(size(G));
            Grest = zeros(size(Grest));
            c = f(u);
            d = 0;
        else
            c = lambda * c + mu1;
            G = lambda * G;
            Grest = lambda * Grest;
            d = mu2;
        end
    end

    function [c, G, Grest, d] = evaluatePolyZonotopeNeuronQuad(obj, c, G, Grest, E, ind, ind_, approx)
        % enclose the hyperbolic tangent activation function with a polynomial
        % zonotope by fitting a quadratic function

        % hyperbolic tangent activation function
        f = obj.f;
        df = obj.df;

        % compute lower and upper bound
        [l, u] = nnHelper.compBoundsPolyZono(c, G, Grest, E, ind, ind_, approx);

        % compute quadratic function that best fits the activation function
        x = linspace(l, u, 10);
        y = f(x);
        temp = nnHelper.leastSquarePolyFunc(x, y, 2);
        c_c = temp(1);
        c_b = temp(2);
        c_a = temp(3);

        % compute the difference between activation function and quad. fit
        L = nnHelper.minMaxDiffQuad(c_a, c_b, c_c, l, u, f, df);
        d = rad(L);

        % evaluate the quadratic approximation on the polynomial zonotope
        [c, G, Grest] = nnHelper.quadApproxPolyZono(c, G, Grest, c_a, c_b, c_c);
        c = c + center(L);
    end

    function [c, G, Grest, d] = evaluatePolyZonotopeNeuronCub(obj, c, G, Grest, E, ind, ind_, approx)
        % enclose the hyperbolic tangent activation function with a polynomial
        % zonotope by fitting a cubic function

        % sigmoid activation function
        f = obj.f;
        df = obj.df;

        % compute lower and upper bound
        [l, u] = nnHelper.compBoundsPolyZono(c, G, Grest, E, ind, ind_, approx);

        % compute cubic function that best fits the activation function
        x = linspace(l, u, 10);
        y = f(x);
        temp = nnHelper.leastSquarePolyFunc(x, y, 3);
        c_d = temp(1);
        c_c = temp(2);
        c_b = temp(3);
        c_a = temp(4);

        % compute the difference between activation function and cubic fit
        L = nnHelper.minMaxDiffCub(c_a, c_b, c_c, c_d, l, u, f, df);
        d = rad(L);

        % evaluate the cubic approximation on the polynomial zonotope
        [c, G, Grest] = nnHelper.cubApproxPolyZono(c, G, Grest, c_a, c_b, c_c, c_d);
        c = c + center(L);
    end

    function r = evaluateTaylmNeuron(obj, input, evParams)
        % enclose the hyperbolic tangent activation function with a Taylor model

        % activation function
        f = obj.f;
        df = obj.getDf(1);
        df2 = obj.getDf(2);
        df3 = obj.getDf(3);

        % compute lower and upper bound
        int = interval(input);
        l = infimum(int);
        u = supremum(int);
        c = center(int);

        % consider different polynomial orders
        if strcmp(evParams.polynomial_approx, 'lin')

            % compute linear enclosure according to Theorem 3.2 in [1]
            lambda = min(df(l), df(u));
            mu1 = 0.5 * (f(u) + f(l) - lambda * (u + l));
            mu2 = 0.5 * (f(u) - f(l) - lambda * (u - l));

            r = lambda * input + interval(mu1-mu2, mu1+mu2);

        elseif strcmp(evParams.polynomial_approx, 'quad')

            % obtain quadratic approximation using the Taylor series expansion
            % according to Example 1 in [2]
            c_a = 0.5 * df2(c);
            c_b = df(c) - df2(c) * c;
            c_c = f(c) - df(c) * c + 0.5 * df2(c) * c^2;

            % compute the difference between activation function and quad. fit
            L = nnHelper.minMaxDiffQuad(c_a, c_b, c_c, l, u, f, df);

            % compute resulting Taylor model
            r = c_a * input^2 + c_b * input + c_c + L;

        elseif strcmp(evParams.polynomial_approx, 'cub')

            % compute cubic approximation using the Taylor series expansion
            c_a = 1 / 6 * df3(c);
            c_b = 0.5 * df2(c) - 0.5 * df3(c) * c;
            c_c = df(c) - df2(c) * c + 0.5 * c^2 * df3(c);
            c_d = f(c) - df(c) * c + 0.5 * df2(c) * c^2 - 1 / 6 * df3(c) * c^3;

            % compute the difference between activation function and cubic fit
            L = nnHelper.minMaxDiffCub(c_a, c_b, c_c, c_d, l, u, f, df);

            % compute resulting Taylor model
            r = c_a * input^3 + c_b * input^2 + c_c * input + c_d + L;
        end
    end

    function [c, G, C, d, l, u] = evaluateConZonotopeNeuron(obj, c, G, C, d, l, u, j, options, evParams)
        throw(CORAerror('CORA:notSupported'));
    end

    function der1 = getDerBounds(obj, l, u)
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

        % return interval
        der1 = interval(df_l, df_u);
    end
end

end

%------------- END OF CODE --------------