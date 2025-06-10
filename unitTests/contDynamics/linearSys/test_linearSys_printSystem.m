function res = test_linearSys_printSystem
% test_linearSys_printSystem - unit test function of printSystem
%
% Syntax:
%    res = test_linearSys_printSystem
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
A = [-2 0; 1 -3];
B = [1; 1];
C = [1 0];
sys = linearSys(A,B,[],C);

printSystem(sys);
printSystem(sys,'high');
printSystem(sys,'high',true);
printSystem(sys,'high',false);

% test fid
filename = 'test.txt';
printSystem(filename,sys,'high',true);
sys_copy = eval(fileread(filename));
assert(isequal(sys,sys_copy));
delete(filename)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
