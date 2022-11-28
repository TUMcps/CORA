function [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotopeQuad(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
% evaluatePolyZonotopeQuad - evaluates the activation layer on a 
%    polyZonotope using polynomials of quadratic order
%
% Syntax:
%    [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotopeQuad(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
%
% Inputs:
%    c, G, Grest, expMat, id, id_, ind, ind_ - parameters of polyZonotope
%    evParams - parameter for NN evaluation
%
% Outputs:
%    updated [c, G, Grest, expMat, id, id_, ind, ind_]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnLayer

% Author:        Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 05-April-2022 (TL)

%------------- BEGIN CODE --------------

num_neurons = size(G, 1);

% for image inputs the input dimension is very large, so we
% perform an order reduction here to reduce comp. time
if ~isempty(evParams.num_generators) && evParams.i == 2 && ...
        size(G, 2) + size(Grest, 2) > evParams.num_generators
    [c, G, Grest, expMat, id] = nnHelper.initialOrderReduction(c, ...
        G, Grest, expMat, id, id_, evParams.num_generators);
    ind = find(prod(ones(size(expMat))- ...
        mod(expMat, 2), 1) == 1);
    ind_ = setdiff(1:size(expMat, 2), ind);
end

% initialization
h = size(G, 2);
q = size(Grest, 2);
temp = q + q * h + 0.5 * (q^2 + q);
c_ = zeros(num_neurons, 1);
G_ = zeros(num_neurons, 0.5*(h^2 + h)+h);
Grest_ = zeros(num_neurons, temp);
d = zeros(num_neurons, 1);

% loop over all neurons in the current layer
for j = 1:num_neurons
    [c_(j), G_(j, :), Grest_(j, :), d(j)] = ...
        obj.evaluatePolyZonotopeNeuronQuad(c(j), G(j, :), Grest(j, :), expMat, ...
        ind, ind_, evParams.bound_approx);
end

% update properties
c = c_;
G = G_;
Grest = Grest_(:, sum(abs(Grest_), 1) > 0);
expMat = [expMat, nnHelper.squareExpMat(expMat)];
[expMat, G] = removeRedundantExponents(expMat, G);

% order reduction
[c, G, Grest, expMat, id, d] = nnHelper.reducePolyZono(c, G, Grest, ...
    expMat, id, d, evParams.num_generators);
temp = diag(d);
Grest = [Grest, temp(:, d > 0)];

% update indices of all-even exponents (for zonotope encl.)
ind = find(prod(ones(size(expMat))-mod(expMat, 2), 1) == 1);
ind_ = setdiff(1:size(expMat, 2), ind);

%------------- END OF CODE --------------