function [c, G, C, d, l, u] = conversionConZonoStarSet(cZ)
% conversionConZonoStarSet - convert a constrained zonotope to a star set
% zonotope in the given dimension
%
% Syntax:
%    [c, G, C, d, l, u] = nnHelper.conversionConZonoStarSet(cZ)
%
% Inputs:
%    cZ - constrained zonotope
%
% Outputs:
%    [c, G, C, d, l, u]
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

m = size(cZ.G, 2);

if isempty(cZ.A)

    c = cZ.c;
    G = cZ.G;
    C = [eye(m); -eye(m)];
    d = ones(2*m, 1);
    u = ones(m, 1);
    l = -ones(m, 1);

else
    % compute point satisfying all constraints with pseudo inverse
    p_ = pinv(cZ.A) * cZ.b;

    % compute null-space of constraints
    T = null(cZ.A);

    % transform boundary constraints of the factor hypercube
    A = [eye(m); -eye(m)];
    b = ones(2*m, 1);
    C = A * T;
    d = b - A * p_;
    c = cZ.c + cZ.G * p_;
    G = cZ.G * T;
    int = interval(polytope(C, d));
    l = infimum(int);
    u = supremum(int);
end

end

% ------------------------------ END OF CODE ------------------------------
