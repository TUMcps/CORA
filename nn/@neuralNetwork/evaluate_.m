function r = evaluate_(obj, input, options, idxLayer)
% evaluate_ - compute the output of a neural network for the given input
%     internal use to speed up computation, use neuralNetwork/evaluate
%
% Syntax:
%    res = evaluate_(obj, input, options)
%    res = evaluate_(obj, input, options, idxlayer)
%
% Inputs:
%    obj - object of class neuralNetwork
%    input - input represented as a numeric or set
%    options - options for neural network evaluation 
%    idxLayer - indices of layers that should be evaluated
%
% Outputs:
%    res - output of the neural network
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate

% Authors:       Tobias Ladner
% Written:       21-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 4
    % default: all layers
    idxLayer = 1:length(obj.layers);
end

% evaluate ----------------------------------------------------------------

if isnumeric(input)                                           % numeric ---
    r = aux_evaluateNumeric(obj, input, options, idxLayer);

elseif isa(input, 'interval')                                % interval ---
    r = aux_evaluateInterval(obj, input, options, idxLayer);

                                                % zonotope/polyZonotope ---
elseif isa(input, 'zonotope') || isa(input, 'polyZonotope')
    r = aux_evaluatePolyZonotope(obj, input, options, idxLayer);

elseif isa(input, 'taylm')                                      % taylm ---
    r = aux_evaluateTaylm(obj, input, options, idxLayer);

elseif isa(input, 'conZonotope')                          % conZonotope ---
    r = aux_evaluateConZonotope(obj, input, options, idxLayer);

else                                                            % other ---
    throw(CORAerror('CORA:notSupported', ...
        ['Set representation ', class(input), ' is not supported.']));
end

end


% Auxiliary functions -----------------------------------------------------

function r = aux_evaluateNumeric(obj, input, options, idxLayer)
    % evaluate numeric
    
    r = input;
    for k = idxLayer
        options.nn.layer_k = k;
        layer_k = obj.layers{k};
        % Store input for backpropgation
        if options.nn.train.backprop
            layer_k.backprop.store.input = r;
        end
        r = layer_k.evaluateNumeric(r, options);
    end

end

function r = aux_evaluateInterval(obj, input, options, idxLayer)
    % evaluate interval
    
    r = input;
    for k = idxLayer
        layer_k = obj.layers{k};
        % Store input for backpropgation
        if options.nn.train.backprop
            layer_k.backprop.store.input = r;
        end
        r = layer_k.evaluateInterval(r, options);
    end

end

function r = aux_evaluatePolyZonotope(obj, input, options, idxLayer)
    % evaluate zonotope/polyZonotope
    
    % we only use polyZonotopes internally
    isZonotope = isa(input, 'zonotope');
    if isZonotope
        % transform to polyZonotope
        % and only use independent generators
        input = polyZonotope(input.c, [], input.G, []);
        options.nn.add_approx_error_to_GI = true;
        options.nn.remove_GI = false;
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
        if isempty(id_)
            id_ = 0;
        end

        ind = find(prod(ones(size(E))-mod(E, 2), 1) == 1);
        ind_ = setdiff(1:size(E, 2), ind);

        if options.nn.order_reduction_sensitivity
            % set sensitivity in each layer (used for order reduction)
            obj.calcSensitivity(c);
        end
    
        % iterate over all layers
        for k = idxLayer
            options.nn.layer_k = k;
            layer_k = obj.layers{k};
            [c, G, GI, E, id, id_, ind, ind_] = ...
                layer_k.evaluatePolyZonotope(c, G, GI, E, id, id_, ind, ind_, options);
            obj.propagateBounds(k, options);
        end
    
        r = polyZonotope(c, G, GI, E, id);

    catch ME
        if strcmp(ME.identifier, 'MATLAB:array:SizeLimitExceeded')
            % out of memory
            throw(CORAerror('CORA:outOfMemory', ...
                sprintf(['While processing layer %i. ', ...
                'Try to set options.nn.num_generators to a lower value.'], k), ...
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

function r = aux_evaluateTaylm(obj, input, options, idxLayer)
    % evaluate taylor model
    
    r = input;
    for k = idxLayer
        options.nn.layer_k = k;
        layer_k = obj.layers{k};
        r = layer_k.evaluateTaylm(r, options);
    end

end

function r = aux_evaluateConZonotope(obj, input, options, idxLayer)
    % evaluate constrained zonotope
    
    % convert constrained zonotope to star set
    [c, G, C, d, l, u] = nnHelper.conversionConZonoStarSet(input);
    
    for k = idxLayer
        options.nn.layer_k = k;
        layer_k = obj.layers{k};
        [c, G, C, d, l, u] = ...
            layer_k.evaluateConZonotope(c, G, C, d, l, u, options);
    end
    % convert star set back to constrained zonotope
    r = nnHelper.conversionStarSetConZono(c, G, C, d, l, u);

end

% ------------------------------ END OF CODE ------------------------------
