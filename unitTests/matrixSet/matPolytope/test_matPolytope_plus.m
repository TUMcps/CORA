function res = test_matPolytope_plus
% test_matPolytope_plus - unit test function for plus
% 
% Syntax:
%    res = test_matPolytope_plus
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
% Written:       02-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% empty matrix zonotope
matP = matPolytope();
resvec(end+1,1) = representsa(matP + 1,'emptySet');

% scalar
V = []; V(:,:,1) = 1; V(:,:,2) = -2;
matP = matPolytope(V);
matPplus = matP + 3;
resvec(end+1,1) = compareMatrices(V + 3,matPplus.V);
matPplus = 3 + matP;
resvec(end+1,1) = compareMatrices(V + 3,matPplus.V);

% nx1 vector
V = []; V(:,:,1) = [0; 1; 1]; V(:,:,2) = [1; -1; -2]; V(:,:,3) = [-2; 0; 1];
matP = matPolytope(V);
matPplus = matP + 3;
resvec(end+1,1) = compareMatrices(V + 3,matPplus.V);

% matrix
V = [];
V(:,:,1) = [0 2; 1 -1; 1 -2];
V(:,:,2) = [1 1; -1 0; -2 1];
V(:,:,3) = [-2 0; 0 1; 1 -1];
matP = matPolytope(V);
matPplus = matP + 3;
resvec(end+1,1) = compareMatrices(V + 3,matPplus.V);

% combine results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
