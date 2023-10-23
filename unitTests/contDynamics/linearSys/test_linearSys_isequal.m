function res = test_linearSys_isequal
% test_linearSys_isequal - unit test for equality check
%
% Syntax:
%    res = test_linearSys_isequal
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
% Written:       09-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

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
c_def = [];

% output matrix: q x n
C = [1 1 0 0;
     0 -0.5 0.5 0];

% feedthrough matrix: q x m
D = [0 0 1;
     0 0 0];

% constant input: q x 1
k = [0; 0.02];

% instantiate different linearSys-objects
sys_AB = linearSys(A,B);
sys_ABC = linearSys(A,B,c_def,C);
sys_ABCD = linearSys(A,B,c_def,C,D);
sys_ABcCDk = linearSys(A,B,c,C,D,k);


% compare instantiated systems
if ~isequal(sys_AB,sys_AB)
    res = false;
elseif isequal(sys_AB,sys_ABC)
    res = false;
elseif ~isequal(sys_ABC,sys_ABC)
    res = false;
elseif isequal(sys_AB,sys_ABCD)
    res = false;
elseif isequal(sys_ABCD,sys_ABcCDk)
    res = false;
elseif ~isequal(sys_ABcCDk,sys_ABcCDk)
    res = false;
end

% case that did not work previously...
sys_lin = linearSys([0 1 -1; 1 0 0; 0 1 0],[0;-1;0],[],[0 0 0.05; 0.05 0.05 0]);
f = @(x,u) [x(2)-x(3); x(1)-u(1); x(2)];
g = @(x,u) [0.05*x(3);0.05*x(1)+0.05*x(2)];
sys_nonlin = nonlinearSys('linearSys',f,g);

if ~isequal(sys_lin,sys_nonlin,1e-14)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
