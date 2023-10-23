function refine(obj, varargin)
% refine - refines the network using the maximum error bound
%
% Syntax:
%    res = refine(obj, max_order, type, heuristic, x, verbose, force_bounds, gamma)
%
% Inputs:
%    obj - object of class neuralNetwork
%    max_order - maximum order for refinement
%    type - "all", "naive", "layer"- or "neuron"-wise refinement
%    heuristic - refinement heuristic
%       "approx_error", "sensitivity", "both", "random", "all",
%       'layer_bias"
%    x - point used for sensitivity analysis
%    verbose - whether additional information should be displayed
%    force_bounds - orders at which to re-compute bounds
%    gamma - threshold for neuron-wise refinement
%
% Outputs:
%    -
%
% References:
%    [1] Ladner, T., et al. (2023). Automatic abstraction refinement in
%        neural network verification using sensitivity analysis. HSCC '23:
%        Proceedings of the 26th International Conference on
%        Hybrid Systems: Computation and Control.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate with 'adaptive'

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   14-April-2022 (sensitivity)
%                28-August-2022 (all, random)
%                16-December-2022 (uniformed input check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin > 8
    throw(CORAerror('CORA:tooManyInputArgs', 8));
end

% parse input
[max_order, type, heuristic, x, verbose, force_bounds, gamma] = ...
    setDefaultValues( ...
        {5, 'layer', 'approx_error', [], false, [], 0.1}, ...
        varargin ...
);

% validate parameters
inputArgsCheck({ ...
    {max_order, 'att', {'double'}, {'integer'}}, ...
    {type, 'str', {'all', 'layer', 'neuron', 'naive'}}, ...
    {heuristic, 'str', {'all', 'random', 'approx_error', 'sensitivity', ...
       'both', 'layer_bias'}}, ...
    {x, 'att', {'double'}}, ...
    {verbose, 'att', {'logical'}}, ...
    {force_bounds, 'att', {'double'}, {'integer'}}, ...
    {gamma, 'att', {'double'}}, ...
})

if isempty(x) && (heuristic == "sensitivity" || heuristic == "both")
    throw(CORAerror('CORA:wrongInputInConstructor', ...
        "No point for sensitivity analysis provided."))
end

% get refinable layers
refinable_layers = obj.getRefinableLayers();

% prepare refinement heuristic
if (heuristic == "sensitivity" || heuristic == "both")
    % calculate sensitivity
    obj.calcSensitivity(x);
    for i = 1:length(refinable_layers)
        layer_i = refinable_layers{i};
        if heuristic == "sensitivity"
            layer_i.refine_heu = vecnorm(layer_i.sensitivity, 2, 1)';
        elseif heuristic == "both"
            layer_i.refine_heu = layer_i.refine_heu .* vecnorm(layer_i.sensitivity, 2, 1)';
        end
    end
elseif strcmp(heuristic, "random")
    for i = 1:length(refinable_layers)
        layer_i = refinable_layers{i};
        layer_i.refine_heu = rand(size(layer_i.refine_heu));
    end
elseif strcmp(heuristic, "layer_bias")
    for i = 1:length(refinable_layers)
        layer_i = refinable_layers{i};
        layer_i.refine_heu = layer_i.refine_heu ./ i;
    end
elseif strcmp(heuristic, "all")
    if verbose
        fprintf("Setting type='all', as heuristic was set to 'all'!\n")
    end
    type = "all";
end

% --- ALL REFINEMENT ---
if type == "all" || type == "naive"
    if verbose
        disp("Refining all neurons...")
    end

    count = 0;
    for i = 1:length(refinable_layers)
        layer_i = refinable_layers{i};
        order_i = layer_i.order;
        count = count + sum(order_i+1 <= max_order);
        layer_i.order = min(order_i+1, max_order);
    end

    if verbose
        fprintf("Refined %d neurons!\n", count)
    end

    % --- LAYER-WISE REFINEMENT ---
elseif type == "layer"
    % determine most sensible choice for refinement
    heu = zeros(length(refinable_layers), 3);
    for i = 1:length(refinable_layers)
        layer_i = refinable_layers{i};
        % reduce with norm (TODO try max?)
        heu(i, 1) = norm(layer_i.refine_heu);
        heu(i, 2) = i; % layer i
        heu(i, 3) = max(layer_i.order); % order in layer i
    end

    % filter and sort
    heu = heu(heu(:, 3) < max_order, :);
    heu = heu(heu(:, 1) > 0, :);
    heu = sortrows(heu, 1, "descend");

    if size(heu, 1) > 0
        % refine
        for i = 1:size(heu, 1)
            layer_i = refinable_layers{heu(i, 2)};
            order_i = layer_i.order;
            if order_i < max_order
                if verbose
                    fprintf("Refined layer %d from order %d to %d!\n", heu(i, 2), max(order_i), max(order_i+1))
                end
                order_i = order_i + 1;
                layer_i.order = order_i;

                if any(max(order_i) == force_bounds)
                    % force re-calculation of bounds in all following layers
                    for j = (heu(i, 2) + 1):length(refinable_layers)
                        layer_j = refinable_layers{j};
                        layer_j.l = [];
                        layer_j.u = [];
                    end
                end

                break;

            end
        end
    else
        if verbose
            fprintf("No layers are left to refine! Either max_order=%d reached or not refineable.\n", max_order);
        end
    end

    % --- NEURON-WISE REFINEMENT ---
elseif type == "neuron"
    % determine most sensible choice for refinement
    heu = zeros(0, 4);
    for i = 1:length(refinable_layers)
        layer_i = refinable_layers{i};
        heu_i = layer_i.refine_heu;
        l = size(heu_i, 1);

        heu_i(:, 2) = ones(l, 1) * i; % layer i
        heu_i(:, 3) = (1:l)'; % neuron in layer i
        heu_i(:, 4) = layer_i.order; % order in layer i

        heu = [heu; heu_i];
    end

    % filter and sort
    heu = heu(heu(:, 4) < max_order, :);
    heu = heu(heu(:, 1) > 0, :);
    heu = sortrows(heu, 1, "descend");

    if size(heu, 1) > 0
        M_max = heu(1, 1);
        heu = heu(heu(:, 1) > gamma*M_max, :);

        l_max = heu(1, 2);
        heu = heu(heu(:, 2) == l_max, :);

        % refine
        for i = 1:size(heu, 1)
            layer_i = refinable_layers{heu(i, 2)};
            order_i = layer_i.order(heu(i, 3));
            M = heu(i, 1);

            if verbose
                fprintf("Refined neuron %d from layer %d from order %d to %d!\n", heu(i, 3), heu(i, 2), order_i, order_i+1)
            end

            order_i = order_i + 1;
            layer_i.order(heu(i, 3)) = order_i;

            if any(order_i == force_bounds)
                % force re-calculate bounds in all following layers
                for j = (heu(i, 2) + 1):length(refinable_layers)
                    layer_j = refinable_layers{j};
                    layer_j.l = [];
                    layer_j.u = [];
                end
            end
        end
    else
        if verbose
            fprintf("No neurons are left to refine! Either max_order=%d reached or not refineable.\n", max_order);
        end
    end
else
    error("Unknown refinement type '%s'", type)
end

end

% ------------------------------ END OF CODE ------------------------------
