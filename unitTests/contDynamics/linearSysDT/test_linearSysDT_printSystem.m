function res = test_linearSysDT_printSystem
% test_linearSysDT_printSystem - unit test function of printSystem
%
% Syntax:
%    res = test_linearSysDT_printSystem
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
% See also: none

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test print of a simple system
A = [-0.4 0.6; 0.6 -0.4];
B = [0; 1];
c = [0;0];
C = [1 0];
dt = 0.4;
sys = linearSysDT(A,B,c,C,dt);

printSystem(sys)
printSystem(sys,'high')
printSystem(sys,'high',true)
printSystem(sys,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
