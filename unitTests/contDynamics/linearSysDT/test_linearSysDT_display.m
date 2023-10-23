function res = test_linearSysDT_display
% test_linearSysDT_display - unit test for display functions
%
% Syntax:
%    res = test_linearSysDT_display
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
% Written:       19-November-2021
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
sys_A = linearSys(A,1);
sys_AB = linearSys(A,B);
sys_ABC = linearSys(A,B,[],C);
sys_ABCD = linearSys(A,B,[],C,D);
sys_ABcCDk = linearSys(A,B,c,C,D,k);

% initialize different linearSysDT-objects
dt = 0.05;
sysDT_A = linearSysDT(sys_A,dt)
sysDT_AB = linearSysDT(sys_AB,dt)
sysDT_ABC = linearSysDT(sys_ABC,dt)
sysDT_ABCD = linearSysDT(sys_ABCD,dt)
sysDT_ABcCDk = linearSysDT(sys_ABcCDk,dt)

% code executed correctly
res = true;

% ------------------------------ END OF CODE ------------------------------
