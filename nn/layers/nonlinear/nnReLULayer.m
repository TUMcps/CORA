classdef nnReLULayer < nnLeakyReLULayer
% nnReLULayer - class for ReLU layers
%
% Syntax:
%    obj = nnReLULayer(name)
%
% Inputs:
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Tran, H.-D., et al. "Star-Based Reachability Analysis of Deep
%        Neural Networks", 2019
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner, Sebastian Sigl, Lukas Koller
% Written:       28-March-2022
% Last update:   07-July-2022 (SS, update nnLeakyReLULayer)
%                04-March-2024 (LK, re-implemented evalNumeric & df_i for performance)
% Last revision: 10-August-2022 (renamed)

% ------------------------------ BEGIN CODE -------------------------------

methods
    % constructor
    function obj = nnReLULayer(name)
        if nargin < 1
            name = [];
        end
        % call super class constructor
        obj@nnLeakyReLULayer(0, name)
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})

    % numeric
    function r = evaluateNumeric(obj, input, options)
        r = max(0, input);
    end
    
    % conZonotope
    function [c, G, C, d, l_, u_] = evaluateConZonotopeNeuron(obj, c, G, C, d, l_, u_, j, options)
        % enclose the ReLU activation function with a constrained zonotope based on
        % the results for star sets in [1]

        n = length(c);
        m = size(G, 2);
        M = eye(n);
        M(:, j) = zeros(n, 1);

        % get lower bound
        if options.nn.bound_approx
            c_ = c(j) + 0.5 * G(j, :) * (u_ - l_);
            G_ = 0.5 * G(j, :) * diag(u_-l_);
            l = c_ - sum(abs(G_));
        else
            problem.f = G(j, :);
            problem.Aineq = C;
            problem.bineq = d;
            problem.Aeq = []; problem.beq = [];
            problem.lb = []; problem.ub = [];
            [~, temp] = CORAlinprog(problem);
            l = c(j) + temp;
        end

        % compute output according to Sec. 3.2 in [1]
        if l < 0

            % compute upper bound
            if options.nn.bound_approx
                u = c_ + sum(abs(G_));
            else
                problem.f = -G(j, :);
                problem.Aineq = C;
                problem.bineq = d;
                problem.Aeq = []; problem.beq = [];
                problem.lb = []; problem.ub = [];
                [~, temp] = CORAlinprog(problem);
                u = c(j) - temp;
            end

            if u <= 0
                % l <= u <= 0 -> linear
                c = M * c;
                G = M * G;
            else
                % compute relaxation
                C1 = [zeros(1, m), -1];
                d1 = 0;
                C2 = [G(j, :), -1];
                d2 = -c(j);
                C3 = [-u / (u - l) * G(j, :), 1];
                d3 = u / (u - l) * (c(j) - l);
                C0 = [C, zeros(size(C, 1), 1)];
                d0 = d;
                C = [C0; C1; C2; C3];
                d = [d0; d1; d2; d3];
                temp = zeros(n, 1);
                temp(j) = 1;
                c = M * c;
                G = M * G;
                G = [G, temp];
                l_ = [l_; 0];
                u_ = [u_; u];
            end
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
