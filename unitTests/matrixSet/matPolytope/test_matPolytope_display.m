function res = test_matPolytope_display
% test_matPolytope_display - unit test function for display (only check for
%    runtime errors)
% 
% Syntax:
%    res = test_matPolytope_display
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

res = true;

% empty matrix zonotope
matP = matPolytope()

% scalar
V{1} = 1; V{2} = -2;
matP = matPolytope(V)

% nx1 vector
V{1} = [0; 1; 1]; V{2} = [1; -1; -2]; V{3} = [-2; 0; 1];
matP = matPolytope(V)

% matrix
V{1} = [0 2; 1 -1; 1 -2];
V{2} = [1 1; -1 0; -2 1];
V{3} = [-2 0; 0 1; 1 -1];
matP = matPolytope(V)

% ------------------------------ END OF CODE ------------------------------
