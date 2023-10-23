function example_manual_matrixSet_expm()
% example_manual_matrixSet_expm - example from the manual demonstrating 
% the expm operation of a matrix set as defined in the manual
%
% Syntax:
%   example_manual_matrixSet_expm()
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

% matrix set
C = [0 1;0 -2.5];
D = [0 0;0 0.5];
A = intervalMatrix(C,D);

% matrix exponentiale
A = expm(A)

% ------------------------------ END OF CODE ------------------------------
