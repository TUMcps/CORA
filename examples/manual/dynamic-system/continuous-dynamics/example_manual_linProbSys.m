function example_manual_linProbSys()
% example_manual_linProbSys - example from the manual demonstrating the 
% linProbSys constructor as defined in the manual
%
% Syntax:
%   example_manual_linProbSys()
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
A = [-1 -4; 4 -1];
B = eye(2);
C = 0.7*eye(2);

% linear system
sys = linProbSys('twoDimSys',A,B,C);

% ------------------------------ END OF CODE ------------------------------
