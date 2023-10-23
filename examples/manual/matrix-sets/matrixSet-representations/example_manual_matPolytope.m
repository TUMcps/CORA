function example_manual_matPolytope()
% example_manual_matPolytope - example from the manual demonstrating the 
% matPolytope constructor as defined in the manual
%
% Syntax:
%   example_manual_matPolytope()
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
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% vertices
V{1} = [1 2; 0 1];
V{2} = [1 3; -1 2];

% matrix polytope
mp = matPolytope(V);

% ------------------------------ END OF CODE ------------------------------
