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
%                30-August-2024 (MW, integrate E and F matrices)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
sys = linearSys();


% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
states = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];
inputs = size(B,2);

% constant offset: n x 1
c = 0.05 * [-4; 2; 3; 1];
c_def = zeros(states,1);

% output matrix: q x n
C = [1 1 0 0;
     0 -0.5 0.5 0];
C_def = 1;
outputs = 2;

% feedthrough matrix: q x m
D = [0 0 1;
     0 0 0];
D_def = 0;

% constant input: q x 1
k = [0; 0.02];
k_def_n = zeros(states,1);
k_def_y = zeros(outputs,1);

% disturbance matrix: n x r
E = [1 0.5; 0 -0.5; 1 -1; 0 1];
E_def = zeros(states,1);
stateDist = size(E,2);

% noise matrix: q x s
F = [1; 0.5];
F_def_n = zeros(states,1);
F_def_y = zeros(outputs,1);
outputDist = size(F,2);


% initialize different linearSys-objects
sys_A = linearSys(A);
sys_AB = linearSys(A,B);
sys_ABC = linearSys(A,B,[],C);
sys_ABCD = linearSys(A,B,[],C,D);
sys_ABcCDk = linearSys(A,B,c,C,D,k);
sys_ABcCDkE = linearSys(A,B,c,C,D,k,E);
sys_ABcCDkEF = linearSys(A,B,c,C,D,k,E,F);
list = {sys_A, sys_AB, sys_ABC, sys_ABCD, sys_ABcCDk, sys_ABcCDkE, sys_ABcCDkEF};

% correct solutions
sys_states = [states, states, states, states, states, states, states];
sys_inputs = [1, inputs, inputs, inputs, inputs, inputs, inputs];
sys_outputs = [states, states, outputs, outputs, outputs, outputs, outputs];
sys_stateDist = [1, 1, 1, 1, 1, stateDist, stateDist];
sys_outputDist = [1, 1, 1, 1, 1, 1, outputDist];
is_c_def = [true, true, false, false, false, false, false];
is_C_def = [true, true, false, false, false, false, false];
is_D_def = [true, true, true, false, false, false, false];
is_k_def = [true, true, true, true, false, false, false];
is_E_def = [true, true, true, true, true, false, false];
is_F_def = [true, true, true, true, true, true, false];

% compare properties of instantiated objects with expected values
for j = 1:length(list)
    
    % current system
    sys = list{j};
    
    % check value of system dimension, number of inputs, number of outputs,
    % number of state disturbances, number of output disturbances
    assert(~(sys.nrOfStates ~= sys_states(j) ...
            || sys.nrOfInputs ~= sys_inputs(j) ...
            || sys.nrOfOutputs ~= sys_outputs(j) ...
            || sys.nrOfDisturbances ~= sys_stateDist(j) ...
            || sys.nrOfNoises ~= sys_outputDist(j)),"loop index %d",j);
    
    % check default values for c, C, D, k, E, F
    assert(~((is_c_def(j) && ~all(sys.c == c_def)) ...
            || (is_C_def(j) && ~all(all(sys.C == C_def))) ...
            || (is_D_def(j) && ~all(all(sys.D == D_def))) ...
            || (is_k_def(j) && ((sys_outputs(j) == states && ~all(sys.k == k_def_n)) ...
                                || (sys_outputs(j) == outputs) && ~all(sys.k == k_def_y))) ...
            || (is_E_def(j) && ~all(all(sys.E == E_def))) ...
            || (is_F_def(j) && ((sys_outputs(j) == states && ~all(sys.F == F_def_n)) ...
                                || (sys_outputs(j) == outputs) && ~all(sys.F == F_def_y)))),"loop index %d",j);
end

% all checks ok
res = true;

% ------------------------------ END OF CODE ------------------------------
