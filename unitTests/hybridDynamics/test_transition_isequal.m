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

% Author:       Mark Wetzlinger
% Written:      26-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% assume true
res = true;

% guard set
guard1 = conHyperplane([-1;0],0,[0,1],0);
guard2 = conHyperplane([-2;0],0,[0,1],0);

% reset function
reset1 = struct('A',[1,0;0,-0.75],'c',[0;0]);
reset2 = struct('A',[1,0;0,-0.75],'c',[1;0]);

% target location
target1 = 1;
target2 = 2;

% synchronization label
syncLabel1 = 'A';
syncLabel2 = 'B';


% same transition
if ~isequal(transition(guard1,reset1,target1,syncLabel1),...
            transition(guard1,reset1,target1,syncLabel1))
    res = false;
end

% different guard set
if isequal(transition(guard1,reset1,target1,syncLabel1),...
           transition(guard2,reset1,target1,syncLabel1))
    res = false;
end

% different reset function
if isequal(transition(guard1,reset1,target1,syncLabel1),...
           transition(guard1,reset2,target1,syncLabel1))
    res = false;
end

% different target location
if isequal(transition(guard1,reset1,target1,syncLabel1),...
           transition(guard1,reset1,target2,syncLabel1))
    res = false;
end

% different synchronization label
if isequal(transition(guard1,reset1,target1,syncLabel1),...
           transition(guard1,reset1,target1,syncLabel2))
    res = false;
end

%------------- END OF CODE --------------
