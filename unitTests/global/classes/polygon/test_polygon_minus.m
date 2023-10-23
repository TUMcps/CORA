function res = test_polygon_minus()
% test_polygon_minus - unit test function for polygon/minus
%
% Syntax:
%    res = test_polygon_minus()
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
% Written:       28-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% generate data
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
x = x(ind);
y = y(ind);

% get polygon
pgon = polygon(x,y);

% test scalar
v = 2;
resvec(end+1) = isequal(pgon-v, polygon(x-v,y-v));
resvec(end+1) = isequal(v-pgon, polygon(v-x,v-y));

% test vector
v = [2;-4];
resvec(end+1) = isequal(pgon-v, polygon(x-v(1),y-v(2)));
resvec(end+1) = isequal(v-pgon, polygon(v(1)-x,v(2)-y));

% check empty case
pgon = polygon();
resvec(end+1) = isequal(pgon-2, pgon);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
