classdef (Abstract) nnActivationLayer < nnLayer
% nnActivationLayer - abstract class for non-linear layers
%
% Syntax:
%    obj = nnActivationLayer(name)
%
% Inputs:
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Kochdumper, N., et al. (2022). Open-and closed-loop neural network
%        verification using polynomial zonotopes. 
%        arXiv preprint arXiv:2207.02715.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner, Lukas Koller
% Written:       28-March-2022
% Last update:   01-April-2022 (moved to class folder)
%                16-February-2023 (combined approx_type)
%                03-May-2023 (LK, backprop)
%                30-May-2023 (approx error/output bounds)
%                02-August-2023 (LK, zonotope batch-eval & -backprop)
%                19-August-2023 (zonotope batch-eval: memory optimizations for GPU training)
%                02-February-2024 (LK, better zonotope backpropagation)
% Last revision: 10-August-2022 (renamed)

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = true     % whether the layer is refineable
end

properties
    % function handles

    f                       % function
    df                      % function derivative

    % adaptive refinement

    order = 1               % order of approximation polynomial
    refine_heu              % heuristic for refinement
    do_refinement = true    % whether the layer should be refined

    l = []                  % lower bound of last input
    u = []                  % upper bound of last input

    merged_neurons = []     % network reduction
end

methods
    % constructor
    function obj = nnActivationLayer(name)
        % call super class constructor
        obj@nnLayer(name)

        % init function handles
        obj.f = @(x) obj.evaluateNumeric(x, struct('backprop', false));
        obj.df = obj.getDf(1);
    end
end

% evaluate (element-wise) -------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})

    % sensitivity
    function S = evaluateSensitivity(obj, S, options)
        if options.nn.store_sensitivity_based_neuron_similarity
            % We compute the cosine similarity between the sensitivity 
            % of the neurons; used to avoid splitting similar neurons in
            % neuralNetwork/verify.
            % Normalize the sensitivity of each neuron.
            nS = S./sqrt(sum(S.^2,1));
            % Compute and store the pairwise similarity matrix.
            obj.backprop.store.similarity = ...
                pagemtimes(nS,'transpose',nS,'none');
        end

        % Obtain the stored input.
        x = obj.backprop.store.input;
        % Compute the sensitivity w.r.t. the inputs.
        S = S.*permute(obj.df(x),[3 1 2]);
        
        if options.nn.store_sensitivity
            % Store the gradient (used for the sensitivity computation).
            obj.sensitivity = S;
        end
    end

    % interval
    function bounds = evaluateInterval(obj, bounds, options)
        if options.nn.reuse_bounds
            % save bounds
            if isempty(obj.l) || isempty(obj.u) || ~all(size(bounds) == size(obj.l))
                obj.l = bounds.inf;
                obj.u = bounds.sup;

                % set bounds
            elseif representsa_(bounds,'emptySet',eps) && ...
                    any(isnan(obj.l)) && any(isnan(obj.u))
                bounds = interval(obj.l, obj.u);
            end

            obj.l = max(obj.l, bounds.inf);
            obj.u = min(obj.u, bounds.sup);
        end

        % propagate through layer
        bounds = evaluateInterval@nnLayer(obj, bounds, options);
    end

    % zonotope/polyZonotope
    [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
    [c, G, GI, d] = evaluatePolyZonotopeNeuron(obj, c, G, GI, E, Es, order, ind, ind_, options)
    
    % zonotope batch evaluation
    [c, G] = evaluateZonotopeBatch(obj, c, G, options);

    % taylm
    r = evaluateTaylm(obj, input, options)
    function r = evaluateTaylmNeuron(obj, input, order, options)
        % enclose the ReLU activation function with a Taylor model by
        % fitting a quadratic function

        % compute lower and upper bound
        int = interval(input);
        l = infimum(int);
        u = supremum(int);

        % compute approx poly + error
        [coeffs, d] = computeApproxPoly(obj, l, u, order, options.nn.poly_method);

        % evaluate
        r = coeffs(end) + interval(-d, d);
        for i=1:length(coeffs)-1
            r = r + coeffs(end-i) * input^i;
        end        
    end
    
    % conZonotope
    [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options)
    function [c, G, C, d, l, u] = evaluateConZonotopeNeuron(obj, c, G, C, d, l, u, j, options)
        throw(CORAerror('CORA:nnLayerNotSupported', obj, 'conZonotope'))
    end

    % backprop ------------------------------------------------------------

    function grad_in = backpropNumeric(obj, input, grad_out, options, updateWeights)
        % backpropagte gradient
        grad_in = obj.df(input) .* grad_out;
    end

    function [gl, gu] = backpropIntervalBatch(obj, l, u, gl, gu, options, updateWeights)
        gl = obj.df(l) .* gl;
        gu = obj.df(u) .* gu;
    end

    % zonotope batch evaluation
    [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options, updateWeights);
end

% Auxiliary functions -----------------------------------------------------

methods

    function df_i = getDf(obj, i)
        df_i = nnHelper.lookupDf(obj,i);
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = [];
        nout = [];
    end

    function outputSize = getOutputSize(obj, inputSize)
        outputSize = inputSize; % for most activation functions
    end

    % approximation polynomial + error

    [coeffs, d] = computeApproxPoly(obj, l, u, varargin)

    function [coeffs, d] = computeApproxError(obj, l, u, coeffs)
        % bound approximation error according to [1, Sec. 3.2]

        % compute the difference between activation function and quad. fit
        [df_l,df_u] = obj.getDerBounds(l, u);
        [diffl,diffu] = nnHelper.minMaxDiffOrder(coeffs, l, u, obj.f, df_l,df_u);
        
        % change polynomial s.t. lower and upper error are equal
        diffc = (diffl+diffu)/2;
        coeffs(end) = coeffs(end) + diffc;
        d = diffu-diffc; % error is radius then.
    end

    % find regions with approximating polynomials
    coeffs = findRegionPolys(obj,tol,order,l_max,u_max,pStart,dStart,pEnd,dEnd)

    function fieldStruct = getFieldStruct(obj)
        fieldStruct = struct;
        if ~isempty(obj.merged_neurons)
            fieldStruct.merged_neurons = obj.merged_neurons;
        end
    end
end

methods (Static)
    layer = instantiateFromString(activation)
end

methods (Abstract)
    [df_l,df_u] = getDerBounds(obj, l, u)
end

methods(Access=protected)
    function [coeffs, d] = computeApproxPolyCustom(obj, l, u, order, poly_method)
        % implement custom polynomial computation in subclass
        coeffs = []; d = [];
    end

    function [xs,xs_m] = computeExtremePointsBatch(obj, m, options)
        throw(CORAerror('CORA:nnLayerNotSupported', obj, 'computeExtremePointsBatch'))
    end
end

end

% ------------------------------ END OF CODE ------------------------------
