function completed = example_matPolytope()
% example_matPolytope - example for polytope matrices
%
% Syntax:
%    completed = example_matPolytope()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: matPolytope

% Authors:       Niklas Kochdumper
% Written:       25-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

P1{1} = [1 2; 3 4]; % 1st vertex of matrix polytope P1
P1{2} = [2 2; 3 3]; % 2nd vertex of matrix polytope P1
matP1 = matPolytope(P1); % instantiate matrix polytope P1

P2{1} = [-1 2; 2 -1]; % 1st vertex of matrix polytope P2
P2{2} = [-1 1; 1 -1]; % 2nd vertex of matrix polytope P2
matP2 = matPolytope(P2); % instantiate matrix polytope P2

matP3 = matP1 + matP2 % perform Minkowski addition and display result
matP4 = matP1 * matP2 % compute multiplication and display result

intP = intervalMatrix(matP1) % compute interval matrix and display result

completed = 1;

% ------------------------------ END OF CODE ------------------------------
