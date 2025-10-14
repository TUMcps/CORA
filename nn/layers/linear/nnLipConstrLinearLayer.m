classdef nnLipConstrLinearLayer < nnLinearLayer
% nnLipConstrLinearLayer - class for linear layers where the weights are
%   constrained by a Lipschitz constant.
%
% Syntax:
%    obj = nnLipConstrLinearLayer(W, b, lambda)
%    obj = nnLipConstrLinearLayer(W, b, lambda, name)
%
% Inputs:
%    W - weight matrix
%    b - bias column vector
%    lambda - Lipschitz constant
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Kitouni, O. et al. Expressive monotonic neural networks. (ICLR). 2023
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Lukas Koller
% Written:       18-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    lambda                  % Lipschitz constant
end

methods
    % constructor
    function obj = nnLipConstrLinearLayer(W, varargin)
        % Parse the constructor arguments.
        [b, lambda, name] = setDefaultValues({0, 1, []}, varargin);
        inputArgsCheck({ ...
            {W, 'att', {'numeric', 'gpuArray'}}; ...
            {b, 'att', {'numeric', 'gpuArray'}}; ...
            {lambda, 'att', {'numeric', 'gpuArray'}}; ...
        })

        % Call constructor of super class.
        obj@nnLinearLayer(W,b,name)
        % Set attributes.
        obj.lambda = lambda;
    end
end

% Auxiliary functions -----------------------------------------------------

methods
    function obj = normWeights(obj,options)
        % Normalize the weights with the Lipschitz constant [1,Eq. 10].
        N = computeWeightNorm(obj);
        obj.W = obj.W*N;
        if options.nn.train.backprop
            % Store the weight normalization for the backpropagation.
            obj.backprop.store.W_norm = N;
        end
    end
end

methods (Access=protected)

    function obj = updateGrad(obj, name, grad, options)
        % We have to the constraint the gradient of weights.
        if strcmp(name,'W')
            if isfield(obj.backprop.store,'W_norm')
                % Obtain the stored weight normalization.
                W_norm = obj.backprop.store.W_norm;
            else
                % There is no stored weight normalization.
                W_norm = 1;
            end
            % Backprop through the normalization.
            grad = grad*W_norm';
        end
        % Update the gradient by calling the super class.
        obj.updateGrad@nnLayer(name,grad,options);
        % Normalize the weights.
        obj.normWeights(options);
    end

    function N = computeWeightNorm(obj)
        % Compute normalization of the weights with the Lipschitz 
        % constant [1,Eq. 10].
        N = diag(1./max(1,obj.lambda*sum(abs(obj.W),1)));
    end
end

end

% ------------------------------ END OF CODE ------------------------------
