function res = test_linearSysDT_linearSysDT
% test_linearSysDT_linearSysDT - unit test for constructor, to see if
%    properties are set correctly
%
% Syntax:
%    res = test_linearSysDT_linearSysDT
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

% define linearSysDT ------------------------------------------------------

% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
n = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];
m = size(B,2);

% constant offset: n x 1
c = 0.05 * [-4; 2; 3; 1];
c_def = zeros(n,1);

% output matrix: q x n
C = [1 1 0 0;
     0 -0.5 0.5 0];
C_def = 1;
y = 2;

% feedthrough matrix: q x m
D = [0 0 1;
     0 0 0];
D_def = 0;

% constant input: q x 1
k = [0; 0.02];
k_def_n = zeros(n,1);
k_def_y = zeros(y,1);

% initialize different linearSys-objects
sys_A = linearSys(A,1);
sys_AB = linearSys(A,B);
sys_ABC = linearSys(A,B,[],C);
sys_ABCD = linearSys(A,B,[],C,D);
sys_ABcCDk = linearSys(A,B,c,C,D,k);

% initialize different linearSysDT-objects
dt = 0.05;
sysDT_A = linearSysDT(sys_A,dt);
sysDT_AB = linearSysDT(sys_AB,dt);
sysDT_ABC = linearSysDT(sys_ABC,dt);
sysDT_ABCD = linearSysDT(sys_ABCD,dt);
sysDT_ABcCDk = linearSysDT(sys_ABcCDk,dt);
list = {sysDT_A, sysDT_AB, sysDT_ABC, sysDT_ABCD, sysDT_ABcCDk};


errmsg_c = 'Default value for c wrong';
errmsg_C = 'Default value for C wrong';
errmsg_D = 'Default value for D wrong';
errmsg_k = 'Default value for k wrong';
errmsg_n = 'System dimension wrong';
errmsg_m = 'Input dimension wrong';
errmsg_y = 'Output dimension wrong';

% correct solutions
sys_n = [n, n, n, n, n];
sys_m = [n, m, m, m, m];
sys_y = [n, n, y, y, y];
is_c_def = [true, true, false, false, false];
is_C_def = [true, true, false, false, false];
is_D_def = [true, true, true, false, false];
is_k_def = [true, true, true, true, false];

% compare properties of instantiated objects with expected values
for j = 1:length(list)
    
    % current system
    sys = list{j};
    
    % check value of system dimension
    assertLoop(sys.nrOfStates == sys_n(j),errmsg_n,[],j);
    % check value of number of inputs
    assertLoop(sys.nrOfInputs == sys_m(j),errmsg_m,[],j);
    % check value of number of outputs
    assertLoop(sys.nrOfOutputs == sys_y(j),errmsg_y,[],j);
    
    % check default values
    if is_c_def(j)
        assertLoop(all(sys.c == c_def),errmsg_c,[],j);
    end
    if is_C_def(j)
        assertLoop(all(all(sys.C == C_def)),errmsg_C,[],j);
    end
    if is_D_def(j)
        assertLoop(all(all(sys.D == D_def)),errmsg_D,[],j);
    end
    if is_k_def(j)
        if sys_y(j) == n
        	assertLoop(all(sys.k == k_def_n),errmsg_k,[],j);
        elseif sys_y(j) == y
            assertLoop(all(sys.k == k_def_y),errmsg_k,[],j);
        end
    end
end

% all checks ok
res = true;

end

% ------------------------------ END OF CODE ------------------------------
