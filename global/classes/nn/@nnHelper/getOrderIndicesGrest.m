function [Grest_start, Grest_end, Grest_ext_start, Grest_ext_end] = getOrderIndicesGrest(Grest, G, order)
% getOrderIndicesGrest - calculates the start and end indices for 
%    polynomial evaluation
%
% Syntax:
%    [Grest_start, Grest_end, Grest_ext_start, Grest_ext_end] = nnHelper.getOrderIndicesGrest(Grest, G, order)
%
% Inputs:
%    Grest - ind. generator matrix
%    G - generator matrix
%    order - order of polynomial
%
% Outputs:
%    Grest_start - start indices of Grest^i, i \in 1:order
%    Grest_end - end indices of Grest^i, i \in 1:order
%    Grest_ext_start - start indices of extended Grest^i, i \in 1:order
%    Grest_ext_end - end indices of extended Grest^i, i \in 1:order
%       - where 'extended' refers to Grest5_ext = [Grest2_ext, Grest3_ext, Grest5]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: evaluatePolyZonotopeAdaptive

% Author:        Tobias Ladner
% Written:       06-April-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% init
Grest_start = zeros(1, order);
Grest_end = zeros(1, order);
Grest_ext_start = zeros(1, order);
Grest_ext_end = zeros(1, order);

% init linear terms
n = size(Grest, 2);
Grest_start(1) = 1;
Grest_end(1) = n;
Grest_ext_start(1) = 1;
Grest_ext_end(1) = n;

% get lengths of G
[~, ~, G_ext_start, G_ext_end] = nnHelper.getOrderIndicesG(G, order);

for o = 2:order
    o1 = floor(o/2);
    o2 = ceil(o/2);

    o1_ext_len = Grest_ext_end(o1) - Grest_ext_start(o1) + 1;
    o2_ext_len = Grest_ext_end(o2) - Grest_ext_start(o2) + 1;

    G_o1_ext_len = G_ext_end(o1) - G_ext_start(o1) + 1;
    G_o2_ext_len = G_ext_end(o2) - G_ext_start(o2) + 1;

    if o1 == o2
        n = o1_ext_len;
        o_len = 0.5 * n * (n + 1) + G_o1_ext_len * o2_ext_len;
    else
        o_len = o1_ext_len * o2_ext_len + G_o1_ext_len * o2_ext_len + G_o2_ext_len * o1_ext_len;
    end

    Grest_start(o) = Grest_end(o-1) + 1;
    Grest_end(o) = Grest_start(o) + o_len - 1;

    Grest_ext_start(o) = Grest_ext_end(o-1) + 1;
    Grest_ext_end(o) = Grest_ext_start(o) + o1_ext_len + o2_ext_len + o_len - 1;
end
end

%------------- END OF CODE --------------