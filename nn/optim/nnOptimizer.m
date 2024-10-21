classdef (Abstract) nnOptimizer
% nnOptimizer - abstract class for optimizer
%
% Syntax:
%    optim = nnOptimizer()
%
% Inputs:
%    lr - learning rate
%    lambda - weight decay
%    lrDecayIter - iteration where learning rate is decreased
%    lrDecay - learning rate decay factor
%
% Outputs:
%    optim - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner, Lukas Koller
% Written:       01-March-2023
% Last update:   25-July-2023 (LK, deleteGrad)
%                02-August-2023 (LK, print)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    lr,
    timestep,
    lambda, % weight decay
    lrDecayIter, lrDecay
end

methods
    % constructor
    function optim = nnOptimizer(lr, lambda, lrDecayIter, lrDecay)
        inputArgsCheck({ ...
            {lr, 'att', 'numeric', {'scalar', 'nonnegative'}}; ...
            {lambda, 'att', 'numeric', {'scalar', 'nonnegative'}}; ...
            {lrDecayIter, 'att', 'numeric'}; ...
            {lrDecay, 'att', 'numeric', {'scalar', 'nonnegative'}}...
        })
        optim.lr = lr;
        optim.lambda = lambda;
        optim.lrDecayIter = lrDecayIter;
        optim.lrDecay = lrDecay;
        % initialize timestep
        optim.timestep = 0;
    end

    function optim = step(optim, nn, options, varargin)
        [idxLayer] = setDefaultValues({1:length(nn)}, varargin);
        % Increment timestep.
        optim.timestep = optim.timestep + 1;
        % Decrease learning rate.
        if ismember(optim.timestep,optim.lrDecayIter)
            optim.lr = optim.lr * optim.lrDecay;
        end
        
        % Updates all learnable parameters.
        for i=1:length(nn)
            layer_i = nn.layers{i};
            names = layer_i.getLearnableParamNames();
            for j=1:length(names)
                name = names{j};
                if optim.lambda ~= 0
                    % Apply weight decay.
                    layer_i.backprop.grad.(name) = ...
                        layer_i.backprop.grad.(name) + optim.lambda*layer_i.(name);
                end
                optim.updateParam(layer_i, name);
                % Clear gradient.
                layer_i.backprop.grad.(name) = 0;
            end
        end
    end

    function optim = deleteGrad(optim, nn, options, varargin)
        [idxLayer] = setDefaultValues({1:length(nn)}, varargin);
        % reset timestep
        optim.timestep = 0;
        % delete gradients
        for i=idxLayer
            layeri = nn.layers{i};
            % Reset backpropagation storage.
            layeri.backprop.store = struct;
            % Reset gradients.
            names = layeri.getLearnableParamNames();
            for j=1:length(names)
                layeri.backprop.grad.(names{j}) = 0;
            end
        end
    end

    function s = print(optim)
        s = sprintf('Optimizer, Learning Rate: %.2e',optim.lr);
    end
end

methods  (Access=protected, Abstract)
    optim = updateParam(optim, nnLayer, name, options)
end

end

% ------------------------------ END OF CODE ------------------------------
