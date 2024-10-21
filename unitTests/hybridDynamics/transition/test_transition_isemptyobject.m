function res = test_transition_isemptyobject
% test_transition_isemptyobject - unit test for emptiness check
%
% Syntax:
%    res = test_transition_isemptyobject
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
% Written:       15-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty object
assert(isemptyobject(transition()));

% transition
guard = polytope([0,1],0,[1,0],0);
reset = linearReset([0, 0; 0, 0.2]);
target = 2;
syncLabel = 'on';

trans = transition(guard,reset,target,syncLabel);
assert(~isemptyobject(trans));
trans = transition(guard,reset,target);
assert(~isemptyobject(trans));

% array of transitions
assert(all(isemptyobject([transition(),transition(guard,reset,target)]) == [true,false]));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
