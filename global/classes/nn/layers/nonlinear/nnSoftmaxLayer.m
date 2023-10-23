classdef nnSoftmaxLayer < nnActivationLayer
% nnSoftmaxLayer - class for Softmax layer
%
% Syntax:
%    obj = nnSoftmaxLayer(name)
%
% Inputs:
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Practical Course SoSe '22 - Report Martina Hinz
%    [2] I. Goodfellow et al. "Deep learning". MIT press., 2016
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: NeuralNetwork

% Authors:       Martina Hinz, Tobias Ladner
% Written:       17-June-2022
% Last update:   ---
% Last revision: 17-August-2022

% ------------------------------ BEGIN CODE -------------------------------

properties
    expLayer
end

methods
    % constructor
    function obj = nnSoftmaxLayer(name)
        if nargin < 2
            name = [];
        end
        % call super class constructor
        obj@nnActivationLayer(name)

        obj.expLayer = nnExpLayer();
    end

    %get i-th derivative
    function df_i = getDf(obj, i)
            % see https://automata88.medium.com/how-to-implement-the-softmax-derivative-independently-from-any-loss-function-ae6d44363a9d
            function r = deriv(s)
                jm = diag(s);
                for i = 1:length(jm)
                    for j = 1:length(jm)
                        if i == j
                            jm(i, j) = s(i) * (1 - s(i));
                        else
                            jm(i, j) = -s(i) * s(j);
                        end
                    end
                end
                r = jm;
            end
            
        df_i = @(x) deriv(x);
    end

    function der1 = getDerBounds(obj, l, u)
        % df_l and df_u as lower and upper bound for the derivative
        der1 = interval(0, 1);
    end

    function tol = minMaxDiffSoftmax(obj, l, u, coeffs_n, der1, dx)
        %tol = 0.0001;
        der2 = nnHelper.getDerInterval(coeffs_n, l, u);
        der2 = -der2; % '-' as we calculate f(x) - p(x)

        der = supremum(abs(der1-der2));

        tol = dx * der;
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})

    % numeric
    function [r, obj] = evaluateNumeric(obj, input, evParams)
        % avoid numerical issues see [2, Chp. 4]
        input = input - max(input);
        z = exp(input);
        n = sum(exp(input), 1);
        r = z ./ n;
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
        % avoid numerical issues
        c = c - max(c);

        % evaluate polyZonotope on expLayer
        poly_method = evParams.poly_method;
        evParams.poly_method = 'singh'; % formerly known as 'lin'
        [c, G, GI, E, id, id_, ind, ind_] = obj.expLayer.evaluatePolyZonotope(c, G, GI, E, id, id_, ind, ind_, evParams); 
        evParams.poly_method = poly_method;
        num_neurons = size(G, 1);
        order = max(obj.order);

        if order > 1
            throw(CORAerror("CORA:wrongValue", 'nnSoftmaxLayer only supports order 1.'))
        end

        % initialization
        [G_start, G_end, ~, ~] = nnHelper.getOrderIndicesG(G, order);
        [~, GI_end, ~, ~] = nnHelper.getOrderIndicesGI(GI, G, order);

        % calculate exponential matrix
        Es = zeros(size(E, 1), G_end(end));
        Es(:, 1:G_end(1)) = E;

        % init
        c_ = zeros(num_neurons, 1);
        G_ = zeros(num_neurons, G_end(end));
        GI_ = zeros(num_neurons, GI_end(end));
        d = zeros(num_neurons, 1);

        if ~all(size(obj.order) == size(c))
            obj.order = ones(size(c)) .* obj.order;
        end

        % sum dimensions
        c_sum = sum(c, 1);
        G_sum = sum(G, 1);
        GI_sum = sum(GI, 1);

        % loop over all neurons in the current layer
        for j = 1:num_neurons
            evParams.j = j;
            order_j = obj.order(j);
            [c_(j), G_j, GI_j, d(j)] = ...
                obj.evaluatePolyZonotopeNeuronSoftmax(c(j), G(j, :), GI(j, :), Es, order_j, ind, ind_, c_sum, G_sum, GI_sum, evParams);

            G_(j, 1:length(G_j)) = G_j;
            GI_(j, 1:length(GI_j)) = GI_j;
        end

        % update properties
        c = c_;
        G = G_;
        E = Es;

        % add approximation error
        GI = GI_(:, sum(abs(GI_), 1) > 0);
        Gd = diag(d);
        Gd = Gd(:, d > 0);
        if evParams.add_approx_error_to_GI

            GI = [GI, Gd];
        else
            G = [G, Gd];
            E = blkdiag(E, eye(size(Gd, 2)));
            id = [id; 1 + (1:size(Gd, 2))' * id_];
            id_ = max(id);
        end
    end

    % zonotope/polyZonotope neuron
    function [c, G, GI, d] = evaluatePolyZonotopeNeuronSoftmax(obj, c, G, GI, Es, order, ind, ind_, c_sum, G_sum, GI_sum, evParams)

        % store information to directly compact generators
        [G_start, G_end, G_ext_start, G_ext_end] = nnHelper.getOrderIndicesG(G, order);

        % extract original exponential map
        E = Es(:, 1:G_end(1));

        % compute lower and upper bounds for both dimensions
        [l_n, u_n] = nnHelper.compBoundsPolyZono(c, G, GI, E, ind, ind_, evParams.bound_approx);
        [l_sum, u_sum] = nnHelper.compBoundsPolyZono(c_sum, G_sum, GI_sum, E, ind, ind_, evParams.bound_approx);

        % der bounds [0,1]
        dx = 0.0001;
        der1 = interval(0, 1); % dx/(l_sum+dx));
        N = 1000000;

        % compute function that best fits the activation function
        x_n = linspace(l_n, u_n, floor(sqrt(N)));
        x_sum = linspace(l_sum, u_sum, floor(sqrt(N)));

        % combine points of both dimensions
        [X_n, X_sum] = meshgrid(x_n, x_sum);
        X_n = reshape(X_n, [], 1);
        X_sum = reshape(X_sum, [], 1);

        % filter out valuse in X_n which are larger than values in
        % X_sum -> not plausible
        idx = X_n <= X_sum;
        X_n = X_n(idx);
        X_sum = X_sum(idx);

        X = [X_n.^(0:order), X_sum.^(0:order)];

        y = X_n ./ X_sum;

        % compute coefficients
        coeffs = pinv(X) * y; % coefficients [1.dim 2.dim]

        % turn around coefficients for polyVal
        co = fliplr(coeffs');
        coeffs_sum = co(:, 1:order+1);
        coeffs_n = co(:, order+2:(order + 1)*2);

        % compute difference fo both dimensions
        diff = y - polyval(coeffs_n, X_n) - polyval(coeffs_sum, X_sum);

        % compute space between points
        dx = (abs(l_n) + abs(u_n)) / sqrt(N);
        tol_n = obj.minMaxDiffSoftmax(l_n, u_n, coeffs_n, der1, dx);

        % compute space between points
        dx = (abs(l_sum) + abs(u_sum)) / sqrt(N);
        tol_sum = obj.minMaxDiffSoftmax(l_sum, u_sum, coeffs_sum, der1, dx);

        L = interval(min(diff)-tol_n-tol_sum, max(diff)+tol_n+tol_sum);

        % approximation error
        d = rad(L);

        sizeC = size(c);
        c_out = zeros(2, sizeC(2));
        sizeG = size(G);
        G_out = zeros(2, sizeG(2));
        sizeGr = size(GI);
        GI_out = zeros(2, sizeGr(2));
        coeffs = reshape(coeffs, [order + 1, 2])';

        % evaluate the approximation on the polynomial zonotope
        c_out(1, :) = coeffs(1, 1) + coeffs(1, 2) * c;
        G_out(1, :) = coeffs(1, 2) * G;
        if ~isempty(GI)
            GI_out(1, :) = coeffs(1, 2) * GI;
        end

        c_out(2, :) = coeffs(2, 1) + coeffs(2, 2) * c_sum;
        G_out(2, :) = coeffs(2, 2) * G_sum;
        if ~isempty(GI_sum)
            GI_out(2, :) = coeffs(2, 2) * GI_sum;
        end

        c = sum(c_out) + center(L);
        G = sum(G_out);
        GI = sum(GI_out);
    end

end
end

% ------------------------------ END OF CODE ------------------------------
