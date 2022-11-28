function r = evaluate(obj, input, evParams)
% evaluate - compute the output of a neural network for the given input
%
% Syntax:
%    res = evaluate(obj,input)
%    res = evaluate(obj,input, evParams)
%
% Inputs:
%    obj - object of class neuralNetwork
%    input - input represented as a numeric or set
%    evParams - struct holding parameters for evaluation
%       .bound_approx: bool whether bounds should be overapproximated
%       .polynomial_approx: how non-linear layers should be approximated
%           - 'lin', 'quad', 'cub'
%       .num_generators: max number of generators for order reduction
%       .add_approx_error_to_Grest: whether 'd' should be added to Grest
%
% Outputs:
%    res - output of the neural network
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Author:       Tobias Ladner
% Written:      28-March-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% validate parameters
if nargin < 3
    evParams = struct;
end

% evParams
if ~isfield(evParams, "bound_approx")
    evParams.bound_approx = false;
end

if ~isfield(evParams, "polynomial_approx")
    evParams.polynomial_approx = "lin";
elseif ~ismember(evParams.polynomial_approx, {'lin', 'quad', 'cub'})
    throw(CORAerror('CORA:wrongFieldValue','evParams.polynomial_approx',...
        {'lin', 'quad', 'cub'}));
end

if ~isfield(evParams, "num_generators")
    evParams.num_generators = [];
end

% execute
if isnumeric(input)
    r = input;
    % evaluate numeric input
    for i = 1:length(obj.layers)
        evParams.i = i;
        layer_i = obj.layers{i};
        r = layer_i.evaluateNumeric(r);
    end

elseif isa(input, 'zonotope')
    % evaluate zonotope input
    Z = input.Z;
    for i = 1:length(obj.layers)
        evParams.i = i;
        layer_i = obj.layers{i};
        Z = layer_i.evaluateZonotope(Z, evParams);
    end
    r = zonotope(Z);

elseif isa(input, 'polyZonotope')
    % evaluate polyZonotope input
    % init values
    c = input.c;
    G = input.G;
    Grest = input.Grest;
    expMat = input.expMat;
    id = input.id;
    id_ = max(id);
    if isempty(G)
        G = zeros(size(c, 1), 0);
    end
    if isempty(Grest)
        Grest = zeros(size(c, 1), 0);
    end
    ind = find(prod(ones(size(input.expMat))-mod(input.expMat, 2), 1) == 1);
    ind_ = setdiff(1:size(input.expMat, 2), ind);

    for i = 1:length(obj.layers)
        evParams.i = i;
        layer_i = obj.layers{i};
        [c, G, Grest, expMat, id, id_, ind, ind_] = ...
            layer_i.evaluatePolyZonotope(c, G, Grest, expMat, id, id_, ind, ind_, evParams);
    end

    r = polyZonotope(c, G, Grest, expMat, id);

elseif isa(input, 'taylm')
    % evaluate taylm input
    r = input;
    for i = 1:length(obj.layers)
        evParams.i = i;
        layer_i = obj.layers{i};
        r = layer_i.evaluateTaylm(r, evParams);
    end

elseif isa(input, 'conZonotope')
    % evaluate conZonotope input

    % convert constrained zonotope to star set
    [c, G, C, d, l, u] = nnHelper.conversionConZonoStarSet(input);

    % predefine options for linear programming for speed-up
    options = optimoptions('linprog', 'display', 'off');

    for i = 1:length(obj.layers)
        evParams.i = i;
        layer_i = obj.layers{i};
        [c, G, C, d, l, u] = ...
            layer_i.evaluateConZonotope(c, G, C, d, l, u, options, evParams);
    end
    % convert star set to constrained zonotope
    r = nnHelper.conversionStarSetConZono(c, G, C, d, l, u);

else
    throw(CORAerror('CORA:notSupported',...
        ['Set representation ' class(input) ' is not supported.']));
end

%------------- END OF CODE --------------