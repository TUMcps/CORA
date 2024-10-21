function res = test_linParamSys_printSystem
% test_linParamSys_printSystem - unit test function of printSystem
%
% Syntax:
%    res = test_linParamSys_printSystem
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
Ac = [-2 0; 1.5 -3]; Aw = [0 0; 0.5 0];
A = intervalMatrix(Ac,Aw);
B = [1; 1];
sys = linParamSys(A,B,'varParam');

printSystem(sys)
printSystem(sys,'high')
printSystem(sys,'high',true)
printSystem(sys,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
