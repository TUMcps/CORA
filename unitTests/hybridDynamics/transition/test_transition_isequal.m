function res = test_transition_isequal
% test_transition_isequal - test function for isequal
%
% Syntax:
%    res = test_transition_isequal
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

% guard set
guard1 = polytope([0,1],0,[-1,0],0);
guard2 = polytope([0,1],0,[-1,1],0);

% reset function
reset1 = linearReset([1,0;0,-0.75],[0;1],[0;0]);
reset2 = linearReset([1,0;0,-0.75],[0;1],[1;0]);

% target location
target1 = 1;
target2 = 2;

% synchronization label
syncLabel1 = 'A';
syncLabel2 = 'B';

% two empty transitions
assert(isequal(transition(),transition()));

% same transition
assert(isequal(transition(guard1,reset1,target1,syncLabel1),...
    transition(guard1,reset1,target1,syncLabel1)));

% same transitions in array
log12 = isequal([transition(guard1,reset1,target1,syncLabel1),...
    transition(guard2,reset2,target2,syncLabel2)],...
    [transition(guard1,reset1,target1,syncLabel1),...
    transition(guard2,reset2,target2,syncLabel2)]);
assert(all(size(log12) == [1,2]));
assert(all(log12));

% different guard set
assert(~isequal(transition(guard1,reset1,target1,syncLabel1),...
    transition(guard2,reset1,target1,syncLabel1)));

% different reset function
assert(~isequal(transition(guard1,reset1,target1,syncLabel1),...
    transition(guard1,reset2,target1,syncLabel1)));

% different target location
assert(~isequal(transition(guard1,reset1,target1,syncLabel1),...
    transition(guard1,reset1,target2,syncLabel1)));

% different synchronization label
assert(~isequal(transition(guard1,reset1,target1,syncLabel1),...
    transition(guard1,reset1,target1,syncLabel2)));

% different array length
assert(~isequal(transition(guard1,reset1,target1,syncLabel1),...
    [transition(guard1,reset1,target1,syncLabel2),...
    transition(guard1,reset1,target1,syncLabel2)]));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
