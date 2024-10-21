function res = test_interval_polygon
% test_interval_polygon - unit test function of polygon
%
% Syntax:
%    res = test_interval_polygon
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
I = interval.empty(2);
pgon = polygon(I);
assert(representsa(pgon,"emptySet"));

% init cases
c = [0;1];
lb = [1;2];
ub = [3;4];

% center
I = interval(c);
pgon = polygon(I);
[res,p] = representsa(pgon,'point');
assert(res & all(withinTol(p,c)));

% full set
I = interval(lb,ub);
pgon = polygon(I);

% conversion is outer-approximative
xs = [I.randPoint(100) I.randPoint(100,'extreme')];
assert(all(pgon.contains(xs,'exact',1e-8)));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
