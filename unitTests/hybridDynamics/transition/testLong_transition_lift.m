function res = testLong_transition_lift
% testLong_transition_lift - test function for projection onto
%    higher-dimensional subspaces
%
% Syntax:
%    res = testLong_transition_lift
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
% Written:       10-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% different guard sets:
% 1D interval
guard_1D = interval(-1,1);
% 2D conHyperplane
guard_2D = conHyperplane([-1;0],0,[0,1],0);
% 3D polytope
guard_3D = polytope([1 0 0],1);
% 4D level set
syms x y x3 vx vy
eq = -y + sin(x);
guard_4D = levelSet(-eq,[x;y;vx;vy],'==');

% reset functions (input not supported!)
% 1D -> 1D nonlinear
reset_1D.f = @(x,u) sin(x(1));
% 2D -> 2D linear without inputs
reset_2D = struct('A',[1,0;0,-0.75],'c',[0;0]);
% 3D -> 3D nonlinear without inputs
reset_3D.f = @(x,u) [x(1)*x(2); -x(3); -x(1)];
% 4D -> 4D linear
reset_4D = struct('A',[1,0,0,0;0,-0.75,0,0;0,0,1,0;0,0,0,1],'c',[0;0;-2;1]);

% target & sync label
target = 1;
syncLabel = 'synclabel';

% instantiate transitions
trans_1D = transition(guard_1D,reset_1D,target);
trans_2D = transition(guard_2D,reset_2D,target,syncLabel);
trans_3D = transition(guard_3D,reset_3D,target);
trans_4D = transition(guard_4D,reset_4D,target,syncLabel);

% project to higher-dimensional space
target_ = 2;

% 1D to [3]D of 4D, with identity
trans_1D_ = lift(trans_1D,4,3,target_,true);

% 2D to [1,3]D of 3D, with identity
trans_2D_ = lift(trans_2D,3,[1,3],target_,true);

% 3D to [2,3,5]D of 6D, with identity
trans_3D_ = lift(trans_3D,6,[2,3,5],target_,true);

% 4D to [1,2,4,5]D of 5D, with identity
trans_4D_ = lift(trans_4D,5,[1,2,4,5],target_,true);

% true solutions
guard_1D_ = interval([-Inf;-Inf;-1;-Inf],[Inf;Inf;1;Inf]);
reset_1D_.f = @(x,u) [x(1);x(2);sin(x(3));x(4)];
trans_1D_true = transition(guard_1D_,reset_1D_,target_);

guard_2D_ = conHyperplane([-1;0;0],0,[0,0,1],0);
reset_2D_ = struct('A',[1,0,0;0,1,0;0,0,-0.75],'c',[0;0;0]);
trans_2D_true = transition(guard_2D_,reset_2D_,target_,syncLabel);

guard_3D_ = polytope([0 1 0 0 0 0],1);
reset_3D_.f = @(x,u) [x(1); x(2)*x(3); -x(5); x(4); -x(2); x(6)];
trans_3D_true = transition(guard_3D_,reset_3D_,target_);

% syms dummy
guard_4D_ = levelSet(-eq,[x;y;x3;vx;vy],'==');
reset_4D_ = struct('A',[1,0,0,0,0;0,-0.75,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1],...
    'c',[0;0;0;-2;1]);
trans_4D_true = transition(guard_4D_,reset_4D_,target_,syncLabel);

% compare solutions
res(end+1,1) = isequal(trans_1D_,trans_1D_true);
res(end+1,1) = isequal(trans_2D_,trans_2D_true);
res(end+1,1) = isequal(trans_3D_,trans_3D_true);
res(end+1,1) = isequal(trans_4D_,trans_4D_true);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
