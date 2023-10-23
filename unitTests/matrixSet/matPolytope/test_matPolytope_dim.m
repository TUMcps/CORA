function res = test_matPolytope_dim
% test_matPolytope_dim - unit test function for dimension read
% 
% Syntax:
%    res = test_matPolytope_dim
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

% Authors:       Mark Wetzlinger
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty matrix zonotope
matP = matPolytope();
res = all(dim(matP) == [0,0]);

% scalar
V{1} = 1; V{2} = -2;
matP = matPolytope(V);
res(end+1,1) = all(dim(matP) == [1,1]);

% nx1 vector
V{1} = [0; 1; 1]; V{2} = [1; -1; -2]; V{3} = [-2; 0; 1];
matP = matPolytope(V);
res(end+1,1) = all(dim(matP) == [3,1]);
res(end+1,1) = dim(matP,1) == 3;

% matrix
V{1} = [0 2; 1 -1; 1 -2];
V{2} = [1 1; -1 0; -2 1];
V{3} = [-2 0; 0 1; 1 -1];
matP = matPolytope(V);
res(end+1,1) = all(dim(matP) == [3,2]);
res(end+1,1) = dim(matP,1) == 3;
res(end+1,1) = dim(matP,2) == 2;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
