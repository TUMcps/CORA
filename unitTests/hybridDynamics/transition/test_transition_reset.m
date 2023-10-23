function res = test_transition_reset
% test_transition_reset - test function for reset of state vector (jump)
%
% Syntax:
%    res = test_transition_reset
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
% Written:       11-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% guard set and target do not matter here, but required for transition obj
guard_1D = polytope(1,0);
guard_2D = interval([-2;-1],[2;5]);
guard_3D = interval([-2;-5;1],[3;2;4]);
target = 1;
syncLabel = 'on';

% different reset functions
reset_1D_lin_auto = struct('A',2,'c',-1);
reset_2D_lin_inputs = struct('A',[1 0;0 2],'B',[2;-1],'c',[-1;0]);
reset_3D_nonlin_auto.f = @(x,u) [x(1); x(2)*x(3); sin(x(2))];

% instantiate transitions
trans_1D = transition(guard_1D,reset_1D_lin_auto,target);
trans_2D = transition(guard_2D,reset_2D_lin_inputs,target,syncLabel);
trans_3D = transition(guard_3D,reset_3D_nonlin_auto,target);

% reset of single point or points (via cell arrays)
x_1D = 1;
x_ = reset(trans_1D,x_1D);
x_true = 1;
res(end+1,1) = withinTol(x_,x_true);

x_2D = {[1;0],[-1;1]}; u_2D = {-1,1};
x_ = reset(trans_2D,x_2D,u_2D);
x_true = {[-2;1],[0;1]};
res(end+1,1) = all(withinTol(x_{1},x_true{1}));
res(end+1,1) = all(withinTol(x_{2},x_true{2}));

x_3D = [1;0;-1];
x_ = reset(trans_3D,x_3D);
x_true = [1;0;0];
res(end+1,1) = all(withinTol(x_,x_true));


% reset of sets
X_1D = interval(0,1);
X_ = reset(trans_1D,X_1D);
X_true = interval(-1,1);
res(end+1,1) = isequal(X_,X_true);

X_2D = interval([-1;0],[1;1]);
U_2D = interval(-1,1);
X_ = reset(trans_2D,X_2D,U_2D);
X_true = interval([-4;-1],[2;3]);
res(end+1,1) = isequal(X_,X_true);

X_3D = {interval([-1;0;0],[0;pi;2]),interval([0;0;0],[1;0;0])};
X_ = reset(trans_3D,X_3D);
X_true = {interval([-1;0;0],[0;2*pi;1]),interval([0;0;0],[1;0;0])};
res(end+1,1) = contains(X_{1},X_true{1});
res(end+1,1) = contains(X_{2},X_true{2});


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
