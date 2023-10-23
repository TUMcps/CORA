function example_manual_linParamSys()
% example_manual_linParamSys - example from the manual demonstrating the 
% linParamSys constructor as defined in the manual
%
% Syntax:
%   example_manual_linParamSys()
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
Ac = [-2 0; 1.5 -3];
Aw = [0 0; 0.5 0];
A = intervalMatrix(Ac,Aw);
B = [1; 1];

% linear parametric system
sys = linParamSys(A,B,'varParam');

% ------------------------------ END OF CODE ------------------------------
