function res = test_matPolytope_isempty
% test_matPolytope_isempty - unit test function for emptiness check
% 
% Syntax:
%    res = test_matPolytope_isempty
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
res = isempty(matP);

% scalar
V{1} = 1; V{2} = -2;
matP = matPolytope(V);
res(end+1,1) = ~isempty(matP);

% nx1 vector
V{1} = [0; 1; 1]; V{2} = [1; -1; -2]; V{3} = [-2; 0; 1];
matP = matPolytope(V);
res(end+1,1) = ~isempty(matP);

% matrix
V{1} = [0 2; 1 -1; 1 -2];
V{2} = [1 1; -1 0; -2 1];
V{3} = [-2 0; 0 1; 1 -1];
matP = matPolytope(V);
res(end+1,1) = ~isempty(matP);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
