function res = test_transition_display
% test_transition_display - test function for display
%
% Syntax:
%    res = test_transition_display
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
% Written:       08-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% empty transition
transition()

% guard sets
guard_hyp = conHyperplane([-1;0],0,[0,1],0);
syms x y vx vy
eq = -y + sin(x);
guard_ls = levelSet(-eq,[x;y;vx;vy],'==');

% reset functions
reset_nonlin.f = @(x,u) [x(1); ...
         x(2); ... 
         ((1-0.8*cos(x(1))^2)*x(3)+1.8*cos(x(1))*x(4))/(1+cos(x(1))^2); ...
         (1.8*cos(x(1))*x(3)+(-0.8+cos(x(1))^2)*x(4))/(1+cos(x(1))^2)];
reset_lin = struct('A',[1,0;0,-0.75],'c',[0;0]);

% target
target = 1;

% sync label
syncLabel = 'synclabel';

% display transition
trans_lin = transition(guard_hyp,reset_lin,target)
trans_lin = transition(guard_hyp,reset_lin,target,syncLabel)
trans_nonlin = transition(guard_ls,reset_nonlin,target)
trans_nonlin = transition(guard_ls,reset_nonlin,target,syncLabel)

% ------------------------------ END OF CODE ------------------------------
