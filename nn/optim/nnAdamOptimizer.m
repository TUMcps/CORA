classdef nnAdamOptimizer < nnOptimizer
% nnAdamOptimizer - adam optimizer
%
% Syntax:
%    optim = nnAdamOptimizer()
%
% Inputs:
%    lr - learning rate
%    beta1 - gradient decay factor, see [1]
%    beta2 - squared gradient decay factor, see [1]
%    epsilon - denominator offset, see [1]
%    lambda - weight decay
%    amsgrad - bool, whether to use amsgrad extension of Adam
%    lrDecayIter - iteration where learning rate is decreased
%    lrDecay - learning rate decay factor
%
% Outputs:
%    optim - generated object
%
% Reference:
%    [1] Kingma, D. and Ba, J. (2015) Adam: A Method for Stochastic 
%        Optimization. Proceedings of the 3rd International Conference on 
%        Learning Representations (ICLR 2015).
%    [2] https://de.mathworks.com/help/deeplearning/ref/nnet.cnn.trainingoptionsadam.html
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       07-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    beta1, beta2 % exponential decay rates for the moment estimates
    epsilon
    amsgrad % bool, whether to use amsgrad extension of Adam
end

methods
    % constructor
    function optim = nnAdamOptimizer(varargin)
        % parse input
        narginchk(0,9)
        [lr, beta1, beta2, epsilon, lambda, amsgrad, ...
            lrDecayIter, lrDecay, gradThreshold] ...
            = setDefaultValues({0.001, 0.9, 0.999, 1e-8, 0, false, ...
                [], 1, 0}, varargin);
        inputArgsCheck({ ...
            {beta1, 'att', 'numeric', {'scalar', 'nonnegative'}}; ...
            {beta2, 'att', 'numeric', {'scalar', 'nonnegative'}}; ...
            {epsilon, 'att', 'numeric', {'scalar', 'nonnegative'}}; ...
            {lambda, 'att', 'numeric', {'scalar', 'nonnegative'}}; ...
            {amsgrad, 'att', 'logical'}; ...
        })

        optim@nnOptimizer(lr, lambda, lrDecayIter, lrDecay, gradThreshold);
        optim.beta1 = beta1;
        optim.beta2 = beta2;
        optim.epsilon = epsilon;
        optim.amsgrad = amsgrad;
    end

    function s = print(optim)
        s = sprintf(['AdamOptimizer, Learning Rate: %.2e, Beta1: %.2e, '...
            'Beta2: %.2e, Epsilon: %.2e, Lambda: %.2e, L2-Grad norm threshold: %.2e'], ...
                optim.lr,optim.beta1,optim.beta2,optim.epsilon,optim.lambda,optim.gradThreshold);
    end
end

methods (Access=protected)
    function optim = updateParam(optim, layer, name, options)
        % Read gradient.
        grad = layer.backprop.grad.(name);

        % Read moment vectors.
        prev_mt = layer.backprop.mt.(name);
        prev_vt = layer.backprop.vt.(name);

        % Compute moment estimates.
        mt = prev_mt + (1 - optim.beta1).*(grad - prev_mt);
        vt = optim.beta2*prev_vt + (1 - optim.beta2).*grad.^2;
        % Store moment estimates.
        layer.backprop.mt.(name) = mt; 
        layer.backprop.vt.(name) = vt;

        if optim.amsgrad
            % Obtain max moment estimate.
            prev_vtMax = layer.backprop.vtMax.(name);
            vt = max(vt,prev_vtMax);
            % Store max moment estimate.
            layer.backprop.vtMax.(name) = vt;   
        end
        % Integrate bias correction into the learning rate.
        stepSize = optim.lr./(1 - optim.beta1.^optim.timestep);
        sqrtVt = sqrt(vt)./sqrt(1 - optim.beta2.^optim.timestep);

        % update weight
        layer.(name) = layer.(name) - stepSize*(mt./(sqrtVt + optim.epsilon));
    end

    function optim = deleteLayerGrad(optim, layer, name, options)
        % Delete additional gradient information.
        layer.backprop.mt.(name) = 0;
        layer.backprop.vt.(name) = 0;
        if optim.amsgrad
            layer.backprop.vtMax.(name) = 0;
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
