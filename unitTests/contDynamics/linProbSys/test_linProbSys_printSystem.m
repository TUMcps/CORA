function res = test_linProbSys_printSystem
% test_linProbSys_printSystem - unit test function of printSystem
%
% Syntax:
%    res = test_linProbSys_printSystem
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
A = [-1 -4; 4 -1];
B = eye(2);
C = 0.7*eye(2);
sys = linProbSys(A,B,C);

printSystem(sys)
printSystem(sys,'high')
printSystem(sys,'high',true)
printSystem(sys,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
