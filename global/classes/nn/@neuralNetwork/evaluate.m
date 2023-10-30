function r = evaluate(obj, input, varargin)
% evaluate - compute the output of a neural network for the given input
%
% Syntax:
%    res = evaluate(obj, input)
%    res = evaluate(obj, input, evParams)
%
% Inputs:
%    obj - object of class neuralNetwork
%    input - input represented as a numeric or set
%    evParams - struct holding parameters for set-based evaluation
%       .bound_approx: bool whether bounds should be overapproximated,
%           or "sample" (not safe!)
%       .reuse_bounds: wheter bounds should be reused
%       .poly_method: how approximation polynomial is found in nonlinear
%           layers, e.g. 'regression', 'singh', 'taylor' ...
%       .num_generators: max number of generators for order reduction
%       .add_approx_error_to_GI: whether 'd' should be added to GI
%       .plot_multi_layer_approx_info: plotting for nnApproximationLayer
%       .max_bounds: max order used in refinement
%       .do_pre_order_reduction: wheter to do a order reduction before
%           evaluating the pZ using the polynomial
%       .max_gens_post: max num_generator before post order reduction
%       .remove_GI: whether to restructue s.t. there remain no ind. gens
%       .force_approx_lin_at: (l, u) distance at which to use 'lin'
%           instead of respective order using 'adaptive'
%       .sort_exponents: whether exponents should be sorted
%       .propagate_bounds: whether bounds should be propagated to next
%           activation layer using interval arithmetic
%       .maxpool_type: for set-based prediction, 'project' or 'regression'
%       .order_reduction_sensitivity: whether sensitivity should be used
%           during order reduction
%       .G: graph object used for graph neural networks
%   idxLayer - indices of layers that should be evaluated
%
% Outputs:
%    res - output of the neural network
%
% References:
%    [1] Kochdumper, N., et al. (2023). Open-and closed-loop neural network
%        verification using polynomial zonotopes. NASA Formal Methods.
%    [2] Ladner, T., et al. (2023). Automatic abstraction refinement in
%        neural network verification using sensitivity analysis. HSCC '23:
%        Proceedings of the 26th International Conference on
%        Hybrid Systems: Computation and Control.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork, nnHelper/validateEvaluateParams

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   29-November-2022 (validateEvaluateParams)
%                16-February-2023 (re-organized structure)
% Last revision: 17-July-2023 (improved readability)

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 2
    throw(CORAerror('CORA:notEnoughInputArgs', 2));
elseif nargin > 5
    throw(CORAerror('CORA:tooManyInputArgs', 5));
end

% validate parameters
[evParams, idxLayer] = setDefaultValues( ...
    {struct, 1:length(obj.layers)}, varargin);

% validate input
inputArgsCheck({ ...
    {obj, 'att', 'neuralNetwork'}; ...
    {input, 'att', {'numeric', 'interval', 'zonotope', 'polyZonotope', ...
        'taylm', 'conZonotope'}}; ...
    {evParams, 'att', 'struct'}; ...
    {idxLayer, 'att', 'numeric', 'vector'}; ...
    })
evParams = nnHelper.validateEvaluateParams(evParams);

% evaluate ----------------------------------------------------------------

if isnumeric(input)                                           % numeric ---
    r = aux_evaluateNumeric(obj, input, evParams, idxLayer);

elseif isa(input, 'interval')                                % interval ---
    r = aux_evaluateInterval(obj, input, evParams, idxLayer);

                                                % zonotope/polyZonotope ---
elseif isa(input, 'zonotope') || isa(input, 'polyZonotope')
    r = aux_evaluatePolyZonotope(obj, input, evParams, idxLayer);

elseif isa(input, 'taylm')                                      % taylm ---
    r = aux_evaluateTaylm(obj, input, evParams, idxLayer);

elseif isa(input, 'conZonotope')                          % conZonotope ---
    r = aux_evaluateConZonotope(obj, input, evParams, idxLayer);

else                                                            % other ---
    throw(CORAerror('CORA:notSupported', ...
        ['Set representation ', class(input), ' is not supported.']));
end

end


% Auxiliary functions -----------------------------------------------------

function r = aux_evaluateNumeric(obj, input, evParams, idxLayer)
    % evaluate numeric
    
    r = input;
    for i = idxLayer
        evParams.i = i;
        layer_i = obj.layers{i};
        r = layer_i.evaluateNumeric(r, evParams);
    end

end

function r = aux_evaluateInterval(obj, input, evParams, idxLayer)
    % evaluate interval
    
    r = input;
    for i = idxLayer
        layer_i = obj.layers{i};
        r = layer_i.evaluateInterval(r, evParams);
    end

end

function r = aux_evaluatePolyZonotope(obj, input, evParams, idxLayer)
    % evaluate zonotope/polyZonotope
    
    % we only use polyZonotopes internally
    isZonotope = isa(input, 'zonotope');
    if isZonotope
        % transform to polyZonotope
        % and only use independent generators
        input = polyZonotope(input.c, [], input.G, []);
        evParams.add_approx_error_to_GI = true;
        evParams.remove_GI = false;
    end
    
    try

        % prepare propagation
        c = input.c;
        G = input.G;
        GI = input.GI;
        E = input.E;
        id = input.id;
        id_ = max(id);
        if isempty(G)
            G = zeros(size(c, 1), 0);
            E = zeros(0,0);
            id = [];
            id_ = 1;
        end
        if isempty(GI)
            GI = zeros(size(c, 1), 0);
        end

        ind = find(prod(ones(size(E))-mod(E, 2), 1) == 1);
        ind_ = setdiff(1:size(E, 2), ind);

        if evParams.order_reduction_sensitivity
            % set sensitivity in each layer (used for order reduction)
            obj.calcSensitivity(c);
        end
    
        % iterate over all layers
        for i = idxLayer
            evParams.i = i;
            layer_i = obj.layers{i};
            [c, G, GI, E, id, id_, ind, ind_] = ...
                layer_i.evaluatePolyZonotope(c, G, GI, E, id, id_, ind, ind_, evParams);
            obj.propagateBounds(i, evParams);
        end
    
        r = polyZonotope(c, G, GI, E, id);

    catch ME
        if strcmp(ME.identifier, 'MATLAB:array:SizeLimitExceeded')
            % out of memory
            throw(CORAerror('CORA:outOfMemory', ...
                sprintf(['While processing layer %i. ', ...
                'Try to set evParams.num_generators to a lower value.'], i), ...
                ME))
        else
            rethrow(ME)
        end
    end

    if isZonotope
        % transform back to zonotope
        r = zonotope(r);
    end

end

function r = aux_evaluateTaylm(obj, input, evParams, idxLayer)
    % evaluate taylor model
    
    r = input;
    for i = idxLayer
        evParams.i = i;
        layer_i = obj.layers{i};
        r = layer_i.evaluateTaylm(r, evParams);
    end

end

function r = aux_evaluateConZonotope(obj, input, evParams, idxLayer)
    % evaluate contrained zonotope
    
    % convert constrained zonotope to star set
    [c, G, C, d, l, u] = nnHelper.conversionConZonoStarSet(input);
    
    % predefine options for linear programming for speed-up
    options = optimoptions('linprog', 'display', 'off');
    
    for i = idxLayer
        evParams.i = i;
        layer_i = obj.layers{i};
        [c, G, C, d, l, u] = ...
            layer_i.evaluateConZonotope(c, G, C, d, l, u, options, evParams);
    end
    % convert star set back to constrained zonotope
    r = nnHelper.conversionStarSetConZono(c, G, C, d, l, u);

end

% ------------------------------ END OF CODE ------------------------------
