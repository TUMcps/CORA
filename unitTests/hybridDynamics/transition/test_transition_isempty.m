function res = test_transition_isempty
% test_transition_isempty - unit test for emptiness check
%
% Syntax:
%    res = test_transition_isempty
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

% Author:       Mark Wetzlinger
% Written:      15-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty object
res = isempty(transition());

% transition
guard = conHyperplane([1,0],0,[0,1],0);
reset = struct('A',[0, 0; 0, 0.2],'c',zeros(2,1));
target = 2;
syncLabel = 'on';

trans = transition(guard,reset,target,syncLabel);
res(end+1,1) = ~isempty(trans);
trans = transition(guard,reset,target);
res(end+1,1) = ~isempty(trans);

% array of transitions
res(end+1,1) = all(isempty([transition(),transition(guard,reset,target)]) == [true,false]);

% combine results
res = all(res);

%------------- END OF CODE --------------