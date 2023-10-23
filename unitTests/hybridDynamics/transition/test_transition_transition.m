function res = test_transition_transition
% test_transition_transition - unit test for constructor of the class
%    transition
%
% Syntax:
%    res = test_transition_transition
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
% Written:       26-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% empty object
trans = transition();

% guard set
guard_2D = conHyperplane([1,0],0,[0,1],0);

% reset function (linear)
reset_2D.A = [0, 0; 0, 0.2]; 
reset_2D.c = zeros(2,1);

% target location
target = 2;

% synchonization label
syncLabel = 'on';

% check standard instantiation
trans = transition(guard_2D,reset_2D,target,syncLabel);

res(end+1,1) = isequal(trans.guard,guard_2D);
res(end+1,1) = all(all(withinTol(trans.reset.A,reset_2D.A))) ...
        && all(withinTol(trans.reset.c,reset_2D.c)) ...
        && trans.reset.stateDim == 2 ...
        && trans.reset.inputDim == 0 ...
        && ~trans.reset.hasInput;
res(end+1,1) = trans.target == target;
res(end+1,1) = strcmp(trans.syncLabel,syncLabel);

% combine results
res = all(res);

% check wrong instantiations
try
    % not enough input arguments
    trans = transition(guard_2D);
    res = false;
end
try
    % not enough input arguments
    trans = transition(guard_2D,reset_2D);
    res = false;
end
try
    % too many input arguments
    trans = transition(guard_2D,reset_2D,target,syncLabel,syncLabel);
    res = false;
end
try
    % guard set is defined using a wrong set representation
    trans = transitions(zonotope([2;1],eye(2)),reset_2D,target);
    res = false;
end
try
    % reset function has wrong fields
    reset_wrong.g = @(x,u) x(1)+u(1);
    trans = transitions(guard_2D,reset_wrong,target);
    res = false;
end
try
    % target not an integer
    trans = transition(guard_2D,reset_2D,1.5);
    res = false;
end
try
    % target negative
    trans = transition(guard_2D,reset_2D,-1);
    res = false;
end
try
    % sync label is numeric
    trans = transition(guard_2D,reset_2D,target,1);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
