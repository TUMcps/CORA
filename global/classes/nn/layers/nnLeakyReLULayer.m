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
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Author:       Sebastian Sigl, Tobias Ladner
% Written:      11-June-2022
% Last update:  ---
% Last revision:10-August-2022 (renamed)

%------------- BEGIN CODE --------------

properties (Constant)
    type = "LeakyReLULayer"
end

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
            name = nnLeakyReLULayer.type;
        end
        % call super class constructor
        obj@nnActivationLayer(name)
        obj.alpha = alpha;
    end

    function df_i = getDf(obj, i)
        if i == 0
            df_i = obj.f;
        elseif i == 1
            df_i = @(x) 1 * (x > 0) + obj.alpha * (x <= 0);
        else
            df_i = @(x) 0;
        end
    end

    function der1 = getDerBounds(obj, l, u)
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
        der1 = interval(df_l, df_u);
    end

    % evaluate
    function r = evaluateNumeric(obj, input)
        r = max(obj.alpha*input, input);
    end

    function [res, d] = evaluateZonotopeNeuron(obj, input)
        % enclose the ReLU activation function with a zonotope according to
        % Theorem 3.1 in [1] here adapted for LeakyReLU

        % compute lower and upper bound
        l = input(1) - sum(abs(input(2:end)));
        u = input(1) + sum(abs(input(2:end)));

        % compute output
        if u < 0
            res = obj.alpha * input;
            d = 0;
        elseif l > 0
            res = input;
            d = 0;
        else
            lambda = (u - obj.alpha * l) / (u - l);
            mu = 0.5 * (u - ((u - obj.alpha * l) / (u - l)) * u);
            res = lambda * input;
            res(1) = res(1) + mu;
            d = mu;
        end
    end

    function [c, G, Grest, d] = evaluatePolyZonotopeNeuronLin(obj, c, G, Grest, E, ind, ind_, approx)
        % enclose the ReLU activation function with a polynomial zonotope according
        % to Theorem 3.1 in [1]

        % compute lower and upper bound
        [l, u] = nnHelper.compBoundsPolyZono(c, G, Grest, E, ind, ind_, approx);

        % compute output
        if u < 0
            % alpha * pZ
            c = obj.alpha * c;
            G = obj.alpha * G;
            Grest = obj.alpha * Grest;
            d = 0;
        elseif l > 0
            % identity
            d = 0;
        else
            lambda = (u - obj.alpha * l) / (u - l);
            mu = 0.5 * (u - ((u - obj.alpha * l) / (u - l)) * u);
            c = lambda * c + mu;
            G = lambda * G;
            Grest = lambda * Grest;
            d = mu;
        end
    end

    function [c, G, Grest, d] = evaluatePolyZonotopeNeuronQuad(obj, c, G, Grest, E, ind, ind_, approx)
        % enclose the LeakyReLU activation function with a polynomial zonotope by
        % fitting a quadratic function

        % properties
        h = length(G);
        q = length(Grest);

        % compute lower and upper bound
        [l, u] = nnHelper.compBoundsPolyZono(c, G, Grest, E, ind, ind_, approx);

        % compute output
        if u < 0 % like ReLU, but scaled
            c = obj.alpha * c;
            G = [G, zeros(1, 0.5*(h^2 + h))];
            Grest = [Grest, zeros(1, q*h+0.5*(q^2 + q))];
            G = obj.alpha * G; % transform identitiy to alpha * identity
            Grest = obj.alpha * Grest;
            d = 0;
        elseif l > 0 % like ReLU
            G = [G, zeros(1, 0.5*(h^2 + h))];
            Grest = [Grest, zeros(1, q*h+0.5*(q^2 + q))];
            d = 0;
        else
            % compute quadratic function that best fits the activation
            % function with f(l)=alpha*l, f'(l)=alpha, f(u)=u
            c_a = u * (1 - obj.alpha) / (u - l)^2;
            c_b = obj.alpha - 2 * c_a * l;
            c_c = c_a * l^2;

            % compute difference between LeakyReLU and quadratic approximation
            L1 = nnHelper.minMaxQuadFun(-c_a, 1-c_b, -c_c, 0, u);
            L2 = nnHelper.minMaxQuadFun(-c_a, obj.alpha-c_b, -c_c, l, 0);
            L = L1 | L2;
            d = rad(L);

            % evaluate the quadratic approximation on the polynomial zonotope
            [c, G, Grest] = nnHelper.quadApproxPolyZono(c, G, Grest, c_a, c_b, c_c);
            c = c + center(L);
        end
    end

    function [c, G, Grest, d] = evaluatePolyZonotopeNeuronCub(obj, c, G, Grest, E, ind, ind_, approx)
        % enclose the LeakyReLU activation function with a polynomial zonotope by
        % fitting a cubic function

        % properties
        h = length(G);
        q = length(Grest);

        % compute lower and upper bound
        [l, u] = nnHelper.compBoundsPolyZono(c, G, Grest, E, ind, ind_, approx);

        % compute output
        if u < 0 % like ReLU, but scaled
            c = obj.alpha * c;
            G = [G, zeros(1, nchoosek(3+h, h)-1-h)];
            temp = nchoosek(3+q, q) - 1 + h * q + 0.5 * (h^2 + h) * q + 0.5 * (q^2 + q) * h - q;
            Grest = [Grest, zeros(1, temp)];
            G = obj.alpha * G;
            Grest = obj.alpha * Grest;
            d = 0;
        elseif l > 0 % like ReLU
            G = [G, zeros(1, nchoosek(3+h, h)-1-h)];
            temp = nchoosek(3+q, q) - 1 + h * q + 0.5 * (h^2 + h) * q + 0.5 * (q^2 + q) * h - q;
            Grest = [Grest, zeros(1, temp)];
            d = 0;
        else
            % compute cubic function that best fits the activation function
            x = linspace(l, u, 10);
            y = max(obj.alpha*x, x);
            temp = nnHelper.leastSquarePolyFunc(x, y, 3);
            c_d = temp(1);
            c_c = temp(2);
            c_b = temp(3);
            c_a = temp(4);

            % compute difference between LeakyReLU and cubic approximation
            L1 = nnHelper.minMaxCubFun(-c_a, -c_b, 1-c_c, -c_d, 0, u);
            L2 = nnHelper.minMaxCubFun(-c_a, -c_b, obj.alpha-c_c, -c_d, l, 0);
            L = L1 | L2;
            d = rad(L);

            % evaluate the cubic approximation on the polynomial zonotope
            [c, G, Grest] = nnHelper.cubApproxPolyZono(c, G, Grest, c_a, c_b, c_c, c_d);
            c = c + center(L);
        end
    end

    function r = evaluateTaylmNeuron(obj, input, evParams)
        % enclose the ReLU activation function with a Taylor model by
        % fitting a quadratic function

        % compute lower and upper bound
        int = interval(input);
        l = infimum(int);
        u = supremum(int);

        % compute output
        if u < 0
            r = obj.alpha * input;
        elseif l > 0
            r = input;
        else

            if strcmp(evParams.polynomial_approx, 'lin')

                % compute linear enclosure according to Theorem 3.1 in [1]
                lambda = (u - obj.alpha * l) / (u - l);
                mu = 0.5 * (u - ((u - obj.alpha * l) / (u - l)) * u);
                r = lambda * input + interval(0, 2*mu);

            elseif strcmp(evParams.polynomial_approx, 'quad')

                % compute quadratic fit to the activation function
                x = linspace(l, u, 10);
                y = max(obj.alpha*x, x);
                temp = nnHelper.leastSquarePolyFunc(x, y, 2);
                c_c = temp(1);
                c_b = temp(2);
                c_a = temp(3);

                % compute difference between ReLU and quadratic approximation
                L1 = nnHelper.minMaxQuadFun(-c_a, 1-c_b, -c_c, 0, u);
                L2 = nnHelper.minMaxQuadFun(-c_a, obj.alpha-c_b, -c_c, l, 0);
                L = L1 | L2;

                % compute resulting Taylor model
                r = c_a * input^2 + c_b * input + c_c + L;

            elseif strcmp(evParams.polynomial_approx, 'cub')

                % compute cubic fit to the activation function
                x = linspace(l, u, 10);
                y = max(obj.alpha*x, x);
                temp = nnHelper.leastSquarePolyFunc(x, y, 3);
                c_d = temp(1);
                c_c = temp(2);
                c_b = temp(3);
                c_a = temp(4);

                % compute difference between ReLU and cubic approximation
                L1 = nnHelper.minMaxCubFun(-c_a, -c_b, 1-c_c, -c_d, 0, u);
                L2 = nnHelper.minMaxCubFun(-c_a, -c_b, obj.alpha-c_c, -c_d, l, 0);
                L = L1 | L2;

                % compute resulting Taylor model
                r = c_a * input^3 + c_b * input^2 + c_c * input + c_d + L;
            end
        end
    end

    function [c, G, C, d, l_, u_] = evaluateConZonotopeNeuron(obj, c, G, C, d, l_, u_, j, options, evParams)
        % TODO LeakyRelu

        throw(CORAerror('CORA:notSupported'));
    end
end

end

%------------- END OF CODE --------------