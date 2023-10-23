function res = test_linearSys_linearSys
% test_linearSys_linearSys - unit test for constructor, to see if
%    properties are set correctly
%
% Syntax:
%    res = test_linearSys_linearSys
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
% Last update:   15-January-2023 (MW, include new syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = false;

% empty case
sys = linearSys();


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
sys_A = linearSys(A);
sys_AB = linearSys(A,B);
sys_ABC = linearSys(A,B,[],C);
sys_ABCD = linearSys(A,B,[],C,D);
sys_ABcCDk = linearSys(A,B,c,C,D,k);
list = {sys_A, sys_AB, sys_ABC, sys_ABCD, sys_ABcCDk};

errmsg_c = 'Default value for c wrong';
errmsg_C = 'Default value for C wrong';
errmsg_D = 'Default value for D wrong';
errmsg_k = 'Default value for k wrong';
errmsg_n = 'System dimension wrong';
errmsg_m = 'Input dimension wrong';
errmsg_y = 'Output dimension wrong';

% correct solutions
sys_n = [n, n, n, n, n];
sys_m = [1, m, m, m, m];
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
    assert(sys.dim == sys_n(j),['Case ' num2str(j) ': ' errmsg_n]);
    % check value of number of inputs
    assert(sys.nrOfInputs == sys_m(j),['Case ' num2str(j) ': ' errmsg_m]);
    % check value of number of outputs
    assert(sys.nrOfOutputs == sys_y(j),['Case ' num2str(j) ': ' errmsg_y]);
    
    % check default values
    if is_c_def(j)
        assert(all(sys.c == c_def),['Case ' num2str(j) ': ' errmsg_c]);
    end
    if is_C_def(j)
        assert(all(all(sys.C == C_def)),['Case ' num2str(j) ': ' errmsg_C]);
    end
    if is_D_def(j)
        assert(all(all(sys.D == D_def)),['Case ' num2str(j) ': ' errmsg_D]);
    end
    if is_k_def(j)
        if sys_y(j) == n
        	assert(all(sys.k == k_def_n),['Case ' num2str(j) ': ' errmsg_k]);
        elseif sys_y(j) == y
            assert(all(sys.k == k_def_y),['Case ' num2str(j) ': ' errmsg_k]);
        end
    end
end

% all checks ok
res = true;

% ------------------------------ END OF CODE ------------------------------
