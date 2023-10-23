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

res = [];

% guard set
guard1 = conHyperplane([-1;0],0,[0,1],0);
guard2 = conHyperplane([-1;1],0,[0,1],0);

% reset function
reset1 = struct('A',[1,0;0,-0.75],'c',[0;0]);
reset2 = struct('A',[1,0;0,-0.75],'c',[1;0]);

% target location
target1 = 1;
target2 = 2;

% synchronization label
syncLabel1 = 'A';
syncLabel2 = 'B';

% two empty transitions
res(end+1,1) = isequal(transition(),transition());

% same transition
res(end+1,1) = isequal(transition(guard1,reset1,target1,syncLabel1),...
    transition(guard1,reset1,target1,syncLabel1));

% same transitions in array
temp = isequal([transition(guard1,reset1,target1,syncLabel1),...
    transition(guard2,reset2,target2,syncLabel2)],...
    [transition(guard1,reset1,target1,syncLabel1),...
    transition(guard2,reset2,target2,syncLabel2)]);
res(end+1,1) = all(size(temp) == [1,2]);
res(end+1,1) = all(temp);

% different guard set
res(end+1,1) = ~isequal(transition(guard1,reset1,target1,syncLabel1),...
    transition(guard2,reset1,target1,syncLabel1));

% different reset function
res(end+1,1) = ~isequal(transition(guard1,reset1,target1,syncLabel1),...
    transition(guard1,reset2,target1,syncLabel1));

% different target location
res(end+1,1) = ~isequal(transition(guard1,reset1,target1,syncLabel1),...
    transition(guard1,reset1,target2,syncLabel1));

% different synchronization label
res(end+1,1) = ~isequal(transition(guard1,reset1,target1,syncLabel1),...
    transition(guard1,reset1,target1,syncLabel2));

% different array length
res(end+1,1) = ~isequal(transition(guard1,reset1,target1,syncLabel1),...
    [transition(guard1,reset1,target1,syncLabel2),...
    transition(guard1,reset1,target1,syncLabel2)]);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
