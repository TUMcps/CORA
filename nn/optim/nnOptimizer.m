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
    gradThreshold % threshold for L2 clipping
end

methods
    % constructor
    function optim = nnOptimizer(lr, lambda, lrDecayIter, lrDecay, gradThreshold)
        inputArgsCheck({ ...
            {lr, 'att', 'numeric', {'scalar', 'nonnegative'}}; ...
            {lambda, 'att', 'numeric', {'scalar', 'nonnegative'}}; ...
            {lrDecayIter, 'att', 'numeric'}; ...
            {lrDecay, 'att', 'numeric', {'scalar', 'nonnegative'}};...
            {gradThreshold, 'att', 'numeric', {'scalar', 'nonnegative'}}...
        })
        optim.lr = lr;
        optim.lambda = lambda;
        optim.lrDecayIter = lrDecayIter;
        optim.lrDecay = lrDecay;
        optim.gradThreshold = gradThreshold;
        % initialize timestep
        optim.timestep = 0;
    end

    function optim = step(optim, nn, options, varargin)
        % Enumerate all layers.
        [layers,~] = nn.enumerateLayers();
        % Parse parameters.
        [idxLayer] = setDefaultValues({1:length(layers)}, varargin);
        % Increment timestep.
        optim.timestep = optim.timestep + 1;
        % Decrease learning rate.
        if ismember(optim.timestep,optim.lrDecayIter)
            optim.lr = optim.lr * optim.lrDecay;
        end
        
        % Updates all learnable parameters.
        for i=idxLayer
            % Obtain the i-th layer.
            layeri = layers{i};
            % Obtain the learnable parameters.
            names = layeri.getLearnableParamNames();
            for j=1:length(names)
                % Obtain the next parameter name.
                name = names{j};
                if optim.lambda ~= 0
                    % Apply weight decay.
                    layeri.backprop.grad.(name) = ...
                        layeri.backprop.grad.(name) + optim.lambda*layeri.(name);
                end
                if optim.gradThreshold > 0
                    % Compute norm of gradient.
                    gradNorm = sqrt(sum(layeri.backprop.grad.(name).^2,'all'));
                    % Scale the gradient.
                    if gradNorm >= optim.gradThreshold
                        gradScale = optim.gradThreshold./gradNorm;
                        layeri.backprop.grad.(name) = ...
                            gradScale.*layeri.backprop.grad.(name);  
                    end
                end
                optim.updateParam(layeri, name);
                % Clear gradient.
                layeri.backprop.grad.(name) = 0;
            end
        end
    end

    function optim = deleteGrad(optim, nn, options, varargin)
        % Enumerate all layers.
        [layers,~] = nn.enumerateLayers();
        % Parse parameters.
        [idxLayer] = setDefaultValues({1:length(layers)}, varargin);
        % Reset timestep.
        optim.timestep = 0;
        % Delete gradients.
        for i=idxLayer
            layeri = layers{i};
            % Reset backpropagation storage.
            layeri.backprop.store = struct;
            % Reset gradients.
            names = layeri.getLearnableParamNames();
            for j=1:length(names)
                % Obtain the j-th name.
                namej = names{j};
                % Delete the stored gradient.
                layeri.backprop.grad.(namej) = 0;
                % Delete additional gradient information, e.g., momentum.
                optim.deleteLayerGrad(layeri,namej);
            end
        end
    end

    function s = print(optim)
        s = sprintf('Optimizer, Learning Rate: %.2e',optim.lr);
    end
end

methods  (Access=protected, Abstract)
    optim = updateParam(optim, layeri, name, options)
    optim = deleteLayerGrad(optim, layeri, name, options)
end

end

% ------------------------------ END OF CODE ------------------------------
