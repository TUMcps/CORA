function res = testLong_transition_synchronizeResets
% testLong_transition_synchronizeResets - test function for synchronization
%    of reset functions due to synchronization labels
%
% Syntax:
%    res = testLong_transition_synchronizeResets
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

% guards, targets
guard1 = conHyperplane([1 1],1,[2 0],0);
guard2 = polytope([1 -1; -1 1],[1;-1]);
target1 = 1;
target2 = 3;

% two linear resets, no inputs
reset1 = struct('A',[2 -1; 1 0],'c',[1;-1]);
reset2 = struct('A',[-3 2; 0 1],'c',[0;-2]);
transList = [transition(guard1,reset1,target1);
            transition(guard2,reset2,target2)];

% global state and input dimensions
n = 2; m = 0;

% synchronize resets
syncReset = synchronizeResets(transList,n,m,[0;0]);

% check result
res = all(all(withinTol(syncReset.A,[-1 1; 1 1])));
res(end+1,1) = all(withinTol(syncReset.c,[1;-3]));


% linear resets with global inputs
reset1 = struct('A',[2 -1; 1 0],'B',[1; 0],'c',[1;-1]);
reset2 = struct('A',[-3 2; 0 1],'c',[0;-2]);
transList = [transition(guard1,reset1,target1);
            transition(guard2,reset2,target2)];

% global state and input dimensions
n = 2; m = 1;

% synchronize resets
syncReset = synchronizeResets(transList,n,m,[0;0]);

% check result
res(end+1,1) = all(all(withinTol(syncReset.A,[-1 1; 1 1])));
res(end+1,1) = all(withinTol(syncReset.B,[1;0]));
res(end+1,1) = all(withinTol(syncReset.c,[1;-3]));


% nonlinear + linear reset, no inputs
reset1 = struct('A',[2 -1; 1 0],'c',[1;-1]);
reset2 = struct('f',@(x,u) [-x(1)*x(2); 0]);
transList = [transition(guard1,reset1,target1);
            transition(guard2,reset2,target2)];

% global state and input dimensions
n = 2; m = 0;

% synchronize resets
syncReset = synchronizeResets(transList,n,m,[0;0]);

% init symbolic variable
x_sym = sym('x',[2;1]);
% plug in synchronized reset
f_val = syncReset.f(x_sym);
% true solution
f_true = @(x) [2*x(1)-x(2)+1 - x(1)*x(2); x(1)-1];
% plug in true solution
f_val_true = f_true(x_sym);

% compare results (diff = 0!)
res(end+1,1) = all(logical(f_val - f_val_true == 0));


% nonlinear resets with global inputs
reset1 = struct('f',@(x,u) [x(2)-u(1); sqrt(x(1))]);
reset2 = struct('f',@(x,u) [-x(1)*x(2); u(1)]);
transList = [transition(guard1,reset1,target1);
            transition(guard2,reset2,target2)];

% global state and input dimensions
n = 2; m = 1;

% synchronize resets
syncReset = synchronizeResets(transList,n,m,[0;0]);

% init symbolic variables
x_sym = sym('x',[2;1]);
u_sym = sym('u',[1;1]);
% plug in synchronized reset
f_val = syncReset.f([x_sym;u_sym]);
% true solution (x3 = u1)
f_true = @(x) [x(2)-x(3) - x(1)*x(2); sqrt(x(1)) + x(3)];
% plug in true solution
f_val_true = f_true([x_sym;u_sym]);

% compare results (diff = 0!)
res(end+1,1) = all(logical(f_val - f_val_true == 0));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
