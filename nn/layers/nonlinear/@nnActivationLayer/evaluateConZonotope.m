function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options)
% evaluateConZonotope - evaluates the activation layer on a conZonotope neuron-wise
%
% Syntax:
%    [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options)
%
% Inputs:
%    c, G, C, d - parameters of conZonotope
%    l, u - bounds
%    options - options for network evaluation (stored in options.nn)
%
% Outputs:
%    updated [c, G, C, d, l, u]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnLayer

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 05-April-2022 (TL)

% ------------------------------ BEGIN CODE -------------------------------

for j = 1:length(c)
    [c, G, C, d, l, u] = obj.evaluateConZonotopeNeuron(c, G, C, d, l, u, j, options);
end

% ------------------------------ END OF CODE ------------------------------
