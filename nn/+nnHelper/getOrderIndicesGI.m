function [GI_start, GI_end, GI_ext_start, GI_ext_end] = getOrderIndicesGI(GI, G, order)
% getOrderIndicesGI - calculates the start and end indices for 
%    polynomial evaluation
%
% Syntax:
%    [GI_start, GI_end, GI_ext_start, GI_ext_end] = nnHelper.getOrderIndicesGI(GI, G, order)
%
% Inputs:
%    GI - ind. generator matrix
%    G - generator matrix
%    order - order of polynomial
%
% Outputs:
%    GI_start - start indices of GI^i, i \in 1:order
%    GI_end - end indices of GI^i, i \in 1:order
%    GI_ext_start - start indices of extended GI^i, i \in 1:order
%    GI_ext_end - end indices of extended GI^i, i \in 1:order
%       - where 'extended' refers to GI5_ext = [GI2_ext, GI3_ext, GI5]
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
GI_start = zeros(1, order);
GI_end = zeros(1, order);
GI_ext_start = zeros(1, order);
GI_ext_end = zeros(1, order);

% init linear terms
n = size(GI, 2);
GI_start(1) = 1;
GI_end(1) = n;
GI_ext_start(1) = 1;
GI_ext_end(1) = n;

% get lengths of G
[~, ~, G_ext_start, G_ext_end] = nnHelper.getOrderIndicesG(G, order);

for o = 2:order
    o1 = floor(o/2);
    o2 = ceil(o/2);

    o1_ext_len = GI_ext_end(o1) - GI_ext_start(o1) + 1;
    o2_ext_len = GI_ext_end(o2) - GI_ext_start(o2) + 1;

    G_o1_ext_len = G_ext_end(o1) - G_ext_start(o1) + 1;
    G_o2_ext_len = G_ext_end(o2) - G_ext_start(o2) + 1;

    if o1 == o2
        n = o1_ext_len;
        o_len = 0.5 * n * (n + 1) + G_o1_ext_len * o2_ext_len;
    else
        o_len = o1_ext_len * o2_ext_len + G_o1_ext_len * o2_ext_len + G_o2_ext_len * o1_ext_len;
    end

    GI_start(o) = GI_end(o-1) + 1;
    GI_end(o) = GI_start(o) + o_len - 1;

    GI_ext_start(o) = GI_ext_end(o-1) + 1;
    GI_ext_end(o) = GI_ext_start(o) + o1_ext_len + o2_ext_len + o_len - 1;
end
end

% ------------------------------ END OF CODE ------------------------------
