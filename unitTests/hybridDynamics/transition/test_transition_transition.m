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

% empty object
trans = transition();

% guard set
guard_2D = polytope([0,1],0,[1,0],0);

% reset function (linear)
reset_2D = linearReset([0, 0; 0, 0.2],[1;0],[0;-1]);

% target location
target = 2;

% synchonization label
syncLabel = 'on';

% check standard instantiation
trans = transition(guard_2D,reset_2D,target,syncLabel);

assert(isequal(trans.guard,guard_2D));
assert(all(withinTol(trans.reset.A,reset_2D.A),"all"));
assert(all(withinTol(trans.reset.c,reset_2D.c)));
assert(trans.reset.preStateDim == 2);
assert(trans.reset.inputDim == 1);
assert(trans.target == target);
assert(strcmp(trans.syncLabel,syncLabel));


% check wrong instantiations
% one input argument, but not copy constructor
assertThrowsAs(@transition,'CORA:wrongValue',guard_2D);

% not enough input arguments
assertThrowsAs(@transition,'CORA:numInputArgsConstructor',guard_2D,reset_2D);

% too many input arguments
assertThrowsAs(@transition,'CORA:numInputArgsConstructor',guard_2D,reset_2D,target,syncLabel,syncLabel);

% guard set is defined using a wrong set representation
assertThrowsAs(@transition,'CORA:wrongValue',zonotope([2;1],eye(2)),reset_2D,target);

% target not an integer
assertThrowsAs(@transition,'CORA:wrongValue',guard_2D,reset_2D,1.5);

% target negative
assertThrowsAs(@transition,'CORA:wrongValue',guard_2D,reset_2D,-1);

% sync label is numeric
assertThrowsAs(@transition,'CORA:wrongValue',guard_2D,reset_2D,target,1);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
