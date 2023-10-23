function example_manual_matrixSet_plus()
% example_manual_matrixSet_plus - example from the manual demonstrating 
% the plus operation of a matrix set as defined in the manual
%
% Syntax:
%   example_manual_matrixSet_plus()
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

% matrix sets
A1 = intervalMatrix([0 1;2 3],[1 2;0 1]);
A2 = intervalMatrix([3 2;2 2],[0 1;1 0]);

% Minkowski addition
res = A1 + A2

% ------------------------------ END OF CODE ------------------------------
