function res = test_linearARX_printSystem
% test_linearARX_printSystem - unit test function of printSystem
%
% Syntax:
%    res = test_linearARX_printSystem
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
A_bar = { [ -0.400 0.600 ; 0.600 -0.400 ] ; [ 0.100 0.000 ; 0.200 -0.500 ] };
B_bar = { [ 0.000 ; 0.000 ] ; [ 0.300 ; -0.700 ] ; [ 0.100 ; 0.000 ] };
dt = 0.100;
sys = linearARX(A_bar,B_bar,dt);

printSystem(sys)
printSystem(sys,'high')
printSystem(sys,'high',true)
printSystem(sys,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
