function [l, u] = compBoundsPolyZono(c, G, GI, E, ind, ind_, approx)
% compBoundsPolyZono - compute the lower and upper bound of a polynomial
% zonotope in the given dimension
%
% Syntax:
%    [l, u] = nnHelper.compBoundsPolyZono(c, G, GI, E, ind, ind_, approx)
%
% Inputs:
%    c - center of polyZonotope in a dimension
%    G - corresponding dep. generator of polyZonotope as row vector
%    GI - corresponding indep. generator of polyZonotope as row vector
%    E - exponential matrix of polyZonotope
%    ind - all even indices
%
% Outputs:
%    [l, u] -  lower and upper bound in the given dimension
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 28-March-2022 (TL)

% ------------------------------ BEGIN CODE -------------------------------

if approx 
    % using zonotope over-approximation
    c_ = c + 0.5 * sum(G(:, ind), 2);

    l = c_ - sum(abs(0.5*G(:, ind)), 2) - sum(abs(G(:, ind_)), 2) - sum(abs(GI), 2);
    u = c_ + sum(abs(0.5*G(:, ind)), 2) + sum(abs(G(:, ind_)), 2) + sum(abs(GI), 2);
else
    % tighter bounds using splitting
    pZ = polyZonotope(c, G, GI, E);
    int = interval(pZ, 'split');
    l = infimum(int);
    u = supremum(int);
end
end

% ------------------------------ END OF CODE ------------------------------
