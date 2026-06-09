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
%    [1] Nolte, N. et al. Expressive monotonic neural networks. (ICLR). 2023
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

methods (Access = {?nnLayer, ?neuralNetwork})

    % Numeric evaluation.
    function r = evaluateNumeric(obj, x, options)
        % Compute the output using the superclass.
        r = evaluateNumeric@nnLinearLayer(obj,x,options);        
        % Compute normalization.
        N = obj.computeWeightNorm();
        % Normalize the output.
        r = N*r;
    end

    % interval 
    function bounds = evaluateInterval(obj, bounds, options)
        % Compute the output using the superclass.
        bounds = evaluateInterval@nnLinearLayer(obj,bounds,options);        
        % Compute normalization.
        N = obj.computeWeightNorm();
        % Normalize the output (N is a diagonal matrix).
        bounds.inf = N*bounds.inf;
        bounds.sup = N*bounds.sup;
    end

    % Zonotope (batch) evaluation.
    function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
        % Compute the output using the superclass.
        [c,G] = evaluateZonotopeBatch@nnLinearLayer(obj,c,G,options);        
        % Compute normalization.
        N = obj.computeWeightNorm();
        % Normalize the output.
        c = N*c;
        G = pagemtimes(N,G);
    end

    % Sensitivity computation.
    function S = evaluateSensitivity(obj, S, options)    
        % Compute normalization.
        N = obj.computeWeightNorm();
        % Backpropagate the sensitivity matrix through the normalization.
        S = pagemtimes(S,N');
        % Compute the gradient using the superclass.
        S = evaluateSensitivity@nnLinearLayer(obj,S,options);
    end

    % Numeric backpropagation.
    function g = backpropNumeric(obj, x, g, options, updateWeights)        
        % Compute normalization.
        N = obj.computeWeightNorm();
        % Backpropagate the gradient through the normalization.
        g = N'*g;
        % Compute the gradient using the superclass.
        g = backpropNumeric@nnLinearLayer(obj,x,g,options,updateWeights);
    end

    % Interval backpropagation.
    function [gl, gu] = backpropIntervalBatch(obj, l, u, gl, gu, ...
            options, updateWeights)
        % Compute normalization.
        N = obj.computeWeightNorm();
        % Backpropagate the gradient through the normalization.
        gl = N'*gl;
        gu = N'*gu;
        % Compute the gradient using the superclass.
        [gl,gu] = backpropIntervalBatch@nnLinearLayer(obj,l,u,gl,gu, ...
            options,updateWeights);
    end

    % Zonotope (batch) backpropagation.
    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, ...
            options, updateWeights)        
        % Compute normalization.
        N = obj.computeWeightNorm();
        % Backpropagate the gradient through the normalization.
        gc = N'*gc;
        gG = pagemtimes(N',gG);
        % Compute the gradient using the superclass.
        [gc,gG] = backpropZonotopeBatch@nnLinearLayer(obj,c,G,gc,gG, ...
            options,updateWeights);
    end
end

% Auxiliary functions -----------------------------------------------------

methods (Access=protected)

    function N = computeWeightNorm(obj)
        % Compute normalization of the weights with the Lipschitz 
        % constant [1,Eq. 10].
        N = diag(1./max(1,obj.lambda^(-1)*sum(abs(obj.W),2)));
        % ... or [1,Eq. 9].
        % N = obj.lambda*diag(1./max(1,sum(abs(obj.W),1)));
    end
end

end

% ------------------------------ END OF CODE ------------------------------
