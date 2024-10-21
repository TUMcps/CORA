function res = test_polygon_mtimes()
% test_polygon_mtimes - unit test function for polygon/mtimes
%
% Syntax:
%    res = test_polygon_mtimes()
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

% get polygon
V = [1 0 0 1; 0 0 1 1];
pgon = polygon(V);

% multiply with matrix
M = [1 2; -3 4];
pgon_res = M * pgon;
V = vertices(pgon_res);
V_true = [ 1 0 2 3 ; -3 0 4 1 ];
assert(compareMatrices(V,V_true))

% left and right multiplication with scalar
M = -2;
pgon_res = M * pgon;
V = vertices(pgon_res);
V_true = [ -2 0 0 -2 ; 0 0 -2 -2 ];
assert(compareMatrices(V,V_true))

M = -2;
pgon_res = pgon * M;
V = vertices(pgon_res);
V_true = [ -2 0 0 -2 ; 0 0 -2 -2 ];
assert(compareMatrices(V,V_true))

% multiplication with higher-dimensional matrices should fail
M = [1 2 3; 1 2 3];
assertThrowsAs(@mtimes,'CORA:notDefined',M,pgon)
assertThrowsAs(@mtimes,'CORA:notDefined',M',pgon)

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
