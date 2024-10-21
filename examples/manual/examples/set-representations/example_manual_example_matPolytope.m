function example_manual_example_matPolytope()
% example_manual_example_matPolytope - example from the manual demonstrating 
%   matrix polytopes
%
% Syntax:
%   example_manual_example_matPolytope()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Tobias Ladner
% Written:       07-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

P1(:,:,1) = [1 2; 3 4]; % 1st vertex of matrix polytope P1
P1(:,:,2) = [2 2; 3 3]; % 2nd vertex of matrix polytope P1
matP1 = matPolytope(P1); % instantiate matrix polytope P1

P2(:,:,1) = [-1 2; 2 -1]; % 1st vertex of matrix polytope P2
P2(:,:,2) = [-1 1; 1 -1]; % 2nd vertex of matrix polytope P2
matP2 = matPolytope(P2); % instantiate matrix polytope P2


matP3 = matP1 + matP2 % perform Minkowski addition and display result
matP4 = matP1*matP2 % compute multiplication of and display result
 
intP = intervalMatrix(matP1) % compute interval matrix and display result

% ------------------------------ END OF CODE ------------------------------
