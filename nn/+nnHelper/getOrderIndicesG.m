function [G_start, G_end, G_ext_start, G_ext_end] = getOrderIndicesG(G, order)
% getOrderIndicesG - calculates the start and end indices for 
%    polynomial evaluation
%
% Syntax:
%    [G_start, G_end, G_ext_start, G_ext_end] = getOrderIndicesG(G, order)
%
% Inputs:
%    G - generator matrix
%    order - order of polynomial
%
% Outputs:
%    G_start - start indices of G^i, i \in 1:order
%    G_end - end indices of G^i, i \in 1:order
%    G_ext_start - start indices of extended G^i, i \in 1:order
%    G_ext_end - end indices of extended G^i, i \in 1:order
%       - where 'extended' refers to G5_ext = [G2_ext, G3_ext, G5]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnActivationLayer/evaluatePolyZonotope

% Authors:       Tobias Ladner
% Written:       06-April-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
G_start = zeros(1, order);
G_end = zeros(1, order);
G_ext_start = zeros(1, order);
G_ext_end = zeros(1, order);

% init linear terms
n = size(G, 2);
G_start(1) = 1;
G_end(1) = n;
G_ext_start(1) = 1;
G_ext_end(1) = n;

for o = 2:order
    o1 = floor(o/2);
    o2 = ceil(o/2);

    o1_ext_len = G_ext_end(o1) - G_ext_start(o1) + 1;
    o2_ext_len = G_ext_end(o2) - G_ext_start(o2) + 1;

    if o1 == o2
        n = o1_ext_len;
        o_len = 0.5 * n * (n + 1);
    else
        o_len = o1_ext_len * o2_ext_len;
    end

    G_start(o) = G_end(o-1) + 1;
    G_end(o) = G_start(o) + o_len - 1;

    G_ext_start(o) = G_ext_end(o-1) + 1;
    G_ext_end(o) = G_ext_start(o) + o1_ext_len + o2_ext_len + o_len - 1;
end
end

% ------------------------------ END OF CODE ------------------------------
