function res = test_linearSys_display
% test_linearSys_display - unit test for display function
%
% Syntax:
%    res = test_linearSys_display
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

% Authors:       Mark Wetzlinger
% Written:       22-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];

% constant offset: n x 1
c = 0.05 * [-4; 2; 3; 1];

% output matrix: q x n
C = [1 1 0 0;
     0 -0.5 0.5 0];

% feedthrough matrix: q x m
D = [0 0 1;
     0 0 0];

% constant input: q x 1
k = [0; 0.02];

% initialize different linearSys-objects
sys_A = linearSys(A,1)
sys_AB = linearSys(A,B)
sys_ABC = linearSys(A,B,[],C)
sys_ABCD = linearSys(A,B,[],C,D)
sys_ABcCDk = linearSys(A,B,c,C,D,k)


% code executed successfully
res = true;

% ------------------------------ END OF CODE ------------------------------
