function r = evaluateZonotope(obj, Z, evParams)
% evaluateZonotope - evaluates the activation layer on a zonotope
%
% Syntax:
%    r = evaluateZonotope(obj, Z, evParams)
%
% Inputs:
%    Z - generator matrix
%    evParams - parameter for NN evaluation
%
% Outputs:
%    r - evaluated zonotope
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

% initialization
num_neurons = size(Z, 1);
Z_ = zeros(num_neurons, size(Z, 2));
d = zeros(num_neurons, 1);

% iterate over all neurons in the current layer
for j = 1:num_neurons
    [Z_(j, :), d(j)] = obj.evaluateZonotopeNeuron(Z(j, :));
end

% order reduction
if ~isempty(evParams.num_generators) && size(Z_, 2) - 1 > evParams.num_generators
    temp = reduce(zonotope(Z_), 'girard', evParams.num_generators/size(Z_, 1));
    Z_ = temp.Z;
end
temp = diag(d);
Z = [Z_, temp(:, d > 0)];
r = Z;

%------------- END OF CODE --------------