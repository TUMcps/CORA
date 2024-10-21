function res = test_polygon_representsa()
% test_polygon_representsa - unit test function for polygon/representsa
%
% Syntax:
%    res = test_polygon_representsa()
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
% See also: polygon

% Authors:       Tobias Ladner
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate data
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
x = x(ind);
y = y(ind);
pgon = polygon(x,y);

% check 'emptySet'
pgon_empty = polygon.empty(2);
assert(representsa(pgon_empty,'emptySet'))
assert(~representsa(pgon,'emptySet'))

% check 'fullspace'
assert(~representsa(pgon,'fullspace'))

% check 'point'
p = [1;2];
pgon_point = polygon(p);
[res,p_res] = representsa(pgon_point,"point");
assert(res & all(withinTol(p,p_res)));
assert(~representsa(pgon,'point'))

% check 'polygon'
assert(representsa(pgon,'polygon'));
assert(representsa(pgon_empty,'polygon'));

% check 'point'
p = [0;0];
pgon_origin = polygon(p);
[res,p_res] = representsa(pgon_origin,"origin");
assert(res & all(withinTol(p,p_res)));
assert(~representsa(pgon,'origin'))
assert(~representsa(pgon_point,'origin'))

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
