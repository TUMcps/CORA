function example_manual_linearSys()
% example_manual_linearSys - example from the manual demonstrating the 
% linearSys constructor as defined in the manual
%
% Syntax:
%   example_manual_linearSys()
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

% Authors:        Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system matrices
A = [-2 0; 1 -3];
B = [1; 1];
C = [1 0];

% linear system
sys = linearSys(A,B,[],C);

% ------------------------------ END OF CODE ------------------------------
