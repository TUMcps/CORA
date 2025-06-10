function res = test_polygon_minkDiff()
% test_polygon_minkDiff - unit test function for polygon/minkDiff
%
% Syntax:
%    res = test_polygon_minkDiff()
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
% Written:       26-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test simple example here, more sophisticated in testLong_polygon_minkDiff

% init pgon
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
pgon = polygon(x(ind),y(ind));

% init subtrahend
I = interval([-1;-1],[1;1])*0.1;

% compute mink diff
pgon_res = minkDiff(pgon,I);

% check point containment
P = [ ...
 0.805, 0.766, 0.713, 0.710, 0.601, 0.504, 0.363, 0.277, 0.568, 0.702, 0.437, 0.481, 0.773, 0.580, 0.804, 0.578, 0.764, 0.738, 0.596, 0.487 ; ...
 0.459, 0.449, 0.193, 0.124, 0.317, 0.138, 0.265, 0.352, 0.535, 0.554, 0.811, 0.770, 0.844, 0.734, 0.456, 0.442, 0.248, 0.132, 0.279, 0.273 ; ...
];
assert(all(contains(pgon_res,P)));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
