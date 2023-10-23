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

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   01-April-2022 (moved to class folder)
%                16-February-2023 (combined approx_type)
%                30-May-2023 (approx error/output bounds)
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
end

methods
    % constructor
    function obj = nnActivationLayer(name)
        % call super class constructor
        obj@nnLayer(name)

        % init function handles
        obj.f = @(x) obj.evaluateNumeric(x);
        obj.df = obj.getDf(1);
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = [];
        nout = [];
    end

    function outputSize = getOutputSize(obj, inputSize)
        outputSize = inputSize; % for most activation functions
    end
end

methods (Static)
    layer = instantiateFromString(activation)
end

methods (Abstract)
    [df_l,df_u] = getDerBounds(obj, l, u)
end

% evaluate (element-wise) -------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})

    % numeric computed in subclass

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, evParams)
       S = S * diag(obj.df(x));
    end

    % interval
    function bounds = evaluateInterval(obj, bounds, evParams)
        if evParams.reuse_bounds
            % save bounds
            if isempty(obj.l) || isempty(obj.u)
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
        bounds = evaluateInterval@nnLayer(obj, bounds, evParams);
    end

    % zonotope/polyZonotope
    [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
    [c, G, GI, d] = evaluatePolyZonotopeNeuron(obj, c, G, GI, E, Es, order, ind, ind_, evParams)
    
    % taylm
    r = evaluateTaylm(obj, input, evParams)
    function r = evaluateTaylmNeuron(obj, input, order, evParams)
        % enclose the ReLU activation function with a Taylor model by
        % fitting a quadratic function

        % compute lower and upper bound
        int = interval(input);
        l = infimum(int);
        u = supremum(int);

        % compute approx poly + error
        [coeffs, d] = computeApproxPoly(obj, l, u, order, evParams.poly_method);

        % evaluate
        r = coeffs(end) + interval(-d, d);
        for i=1:length(coeffs)-1
            r = r + coeffs(end-i) * input^i;
        end        
    end
    
    % conZonotope
    [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options, evParams)
    function [c, G, C, d, l, u] = evaluateConZonotopeNeuron(obj, c, G, C, d, l, u, j, options, evParams)
        throw(CORAerror('CORA:nnLayerNotSupported', obj, 'conZonotope'))
     end
end

% Auxiliary functions -----------------------------------------------------

methods

    function df_i = getDf(obj, i)
        df_i = nnHelper.lookupDf(obj,i);
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
end

methods(Access=protected)
    function [coeffs, d] = computeApproxPolyCustom(obj, l, u, order, poly_method)
        % implement custom polynomial computation in subclass
        coeffs = []; d = [];
    end
end

end

% ------------------------------ END OF CODE ------------------------------
