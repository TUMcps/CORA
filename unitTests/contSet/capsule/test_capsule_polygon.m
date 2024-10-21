function res = test_capsule_polygon
% test_capsule_polygon - unit test function of polygon
%
% Syntax:
%    res = test_capsule_polygon
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
C = capsule.empty(2);
pgon = polygon(C);
assert(representsa(pgon,"emptySet"));

% init cases
c = [1;2];
g = [2;1];
r = 1;

% center
C = capsule(c);
pgon = polygon(C);
[res,p] = representsa(pgon,'point');
assert(res & all(withinTol(p,c)));

% full set
C = capsule(c,g,r);
pgon = polygon(C);

% conversion is inner-approximative
assert(C.contains(pgon));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
