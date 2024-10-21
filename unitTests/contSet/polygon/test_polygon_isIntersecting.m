function res = test_polygon_isIntersecting()
% test_polygon_isIntersecting - unit test function for polygon/isIntersecting
%
% Syntax:
%    res = test_polygon_isIntersecting()
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
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate data
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
x = x(ind);
y = y(ind);

% get polygon
pgon = polygon(x,y);

% check numeric
assert(isIntersecting(pgon,[0.5;0.5]));
assert(isIntersecting(pgon,[x(1);y(1)]));
assert(~isIntersecting(pgon,[2;1]));

% check polygon
assert(isIntersecting(pgon,pgon));

x = [0.1 0.4 0.25 0.0];
y = [0.5 0.6 0.75 0.8];
pgon2 = polygon(x,y);
assert(isIntersecting(pgon,pgon2));

x = [0.1 0.2 0.25 0.0];
y = [0.5 0.6 0.75 0.8];
pgon2 = polygon(x,y);
assert(~isIntersecting(pgon,pgon2));

% check hole
pgon_I = polygon(interval([0.4;0.4],[0.75;0.75]));
pgon_hole = pgon \ pgon_I;

assert(~isIntersecting(pgon_hole,enlarge(pgon_I,0.9)));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
