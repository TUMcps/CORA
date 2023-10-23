function example_manual_linearSysDT()
% example_manual_linearSysDT - example from the manual demonstrating the 
% linearSysDT constructor as defined in the manual
%
% Syntax:
%   example_manual_linearSysDT()
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
A = [-0.4 0.6; 0.6 -0.4];
B = [0; 1];
C = [1 0];

% sampling time
dt = 0.4;

% linear discrete-time system
sys = linearSysDT(A,B,[],C,dt);

% ------------------------------ END OF CODE ------------------------------
