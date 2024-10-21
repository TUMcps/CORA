function cZ = conversionStarSetConZono(c, G, C, d, l_, u_)
% conversionStarSetConZono - convert a star set to a constrained zonotope
%    zonotope in the given dimension
%
% Syntax:
%    cZ = nnHelper.conversionStarSetConZono(c, G, C, d, l_, u_)
%
% Inputs:
%    parameters of StarSet
%
% Outputs:
%    cZ - constraint zonotope
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

% normalize halfspace normal vector length
temp = sqrt(sum(C.^2, 2));
C = diag(temp) * C;
d = temp .* d;

m = size(G, 2);

% find hypercube constraints
u = zeros(m, 1);
ind = find(any(C' == 1));
ind_ = setdiff(1:size(C, 1), ind);
C_ = C(ind, :);
d_ = d(ind);

for i = 1:m
    ind = find(C_(:, i) == 1);
    if isempty(ind)
        u(i) = u_(i);
    else
        u(i) = min(d_(ind));
    end
end

C = C(ind_, :);
d = d(ind_);

l = zeros(m, 1);
ind = find(any(C' == -1));
ind_ = setdiff(1:size(C, 1), ind);
C_ = C(ind, :);
d_ = d(ind);

for i = 1:m
    ind = find(C_(:, i) == -1);
    if isempty(ind)
        l(i) = l_(i);
    else
        l(i) = max(-d_(ind));
    end
end

C = C(ind_, :);
d = d(ind_);

% scale constraints according to hypercube dimensions
int = interval(l, u);
cen = center(int);
R = diag(rad(int));
d = d - C * cen;
C = C * R;
c = c + G * cen;
G = G * R;

% represent inequality constraints as equivalent equality constraints
if ~isempty(C)
    Z = zonotope(interval(-ones(m, 1), ones(m, 1)));
    l = infimum(interval(C*Z));
    int = interval(l, d);
    cen = center(int);
    r = rad(int);
    A = [C, diag(r)];
    b = cen;
    G = [G, zeros(size(G, 1), length(r))];
else
    A = [];
    b = [];
end

% construct resulting constrained zonotope
cZ = conZonotope(c, G, A, b);
end

% ------------------------------ END OF CODE ------------------------------
