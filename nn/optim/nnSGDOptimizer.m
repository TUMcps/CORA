classdef nnSGDOptimizer < nnOptimizer
% nnSGDOptimizer - gradient descent optimizer with optional momentum
%
% Syntax:
%    optim = nnSGDOptimizer()
%
% Inputs:
%    lr - learning rate
%    momentum - momentum
%    lambda - weight decay
%    lrDecayIter - iteration where learning rate is decreased
%    lrDecay - learning rate decay factor
%
% Outputs:
%    optim - generated object
%
% Reference:
%    [1] https://keras.io/api/optimizers/sgd/
%    [2] https://de.mathworks.com/help/deeplearning/ref/trainingoptions.html#bu80qkw-3_head
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner, Lukas Koller
% Written:       01-March-2023
% Last update:   27-April-2023
%                25-July-2023 (LK, implemented deleteGrad)
%                02-August-2023 (LK, added print function)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    momentum
end

methods
    % constructor
    function optim = nnSGDOptimizer(varargin)
        % parse input
        narginchk(0,6)
        [lr, momentum, lambda, lrDecayIter, lrDecay, gradThreshold] = ...
            setDefaultValues({0.001, 0.9, 0, [], 1, 0}, varargin);
        inputArgsCheck({ ...
            {momentum, 'att', 'numeric', {'scalar', 'nonnegative'}}; ...
        })

        optim@nnOptimizer(lr, lambda, lrDecayIter, lrDecay, gradThreshold);
        optim.momentum = momentum;
    end

    function s = print(optim)
        s = sprintf('SGDOptimizer, Learning Rate: %.2e, Momentum: %.2e', ...
            optim.lr,optim.momentum);
    end
end

methods (Access=protected)
    function optim = updateParam(optim, layer, name, options)
        % Read gradient.
        grad = layer.backprop.grad.(name);
        % Read gradient velocity.
        vel = layer.backprop.vel.(name);

        % Update gradient velocity.
        gradUpdate = optim.momentum*vel - optim.lr*grad;
        layer.backprop.vel.(name) = gradUpdate;

        % Update weight.
        layer.(name) = layer.(name) + gradUpdate;
    end

    function optim = deleteLayerGrad(optim, layer, name, options)
        % Delete additional gradient information.
        layer.backprop.vel.(name) = 0;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
