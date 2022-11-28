function [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotopeLin(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
% evaluatePolyZonotopeLin - evaluates the activation layer on a 
%    polyZonotope using polynomials of linear order
%
% Syntax:
%    [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotopeLin(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
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
% initialization
c_ = zeros(num_neurons, 1);
G_ = zeros(num_neurons, size(G, 2));
Grest_ = zeros(num_neurons, size(Grest, 2));
d = zeros(num_neurons, 1);

% loop over all neurons in the current layer
for j = 1:num_neurons
    % calculate approximation in corresponding Subclass
    [c_(j), G_(j, :), Grest_(j, :), d(j)] = ...
        obj.evaluatePolyZonotopeNeuronLin(c(j), G(j, :), Grest(j, :), expMat, ind, ind_, evParams.bound_approx);
end

% update properties
c = c_;
G = G_;
Grest = Grest_;

% order reduction
[c, G, Grest, expMat, id, d] = nnHelper.reducePolyZono(c, G, Grest, ...
    expMat, id, d, evParams.num_generators);
temp = diag(d);
Grest = [Grest, temp(:, d > 0)];

% update indices of all-even exponents (for zonotope encl.)
if size(G, 2) < size(G_, 2)
    ind = find(prod(ones(size(expMat))- ...
        mod(expMat, 2), 1) == 1);
    ind_ = setdiff(1:size(expMat, 2), ind);
end

%------------- END OF CODE --------------