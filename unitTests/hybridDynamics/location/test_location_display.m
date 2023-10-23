function res = test_location_display
% test_location_display - test function for display
%
% Syntax:
%    res = test_location_display
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

% empty location
location()

% dynamics
f = @(x,u) [x(3); x(4); 0; -9.81+0.01*x(4)^2];
sys_nonlin = nonlinearSys(f);
sys_lin = linearSys([-4 1; 1 -4],[],[2;0]);

% invariant
inv_poly = polytope([-1,0],0);
syms x y vx vy
eq = -y + sin(x);
inv_ls = levelSet(eq,[x;y;vx;vy],'<=');

% guard set
guard_hyp = conHyperplane([-1;0],0,[0,1],0);
guard_ls = levelSet(-eq,[x;y;vx;vy],'==');

% reset function
reset_nonlin.f = @(x,u) [x(1); ...
         x(2); ... 
         ((1-0.8*cos(x(1))^2)*x(3)+1.8*cos(x(1))*x(4))/(1+cos(x(1))^2); ...
         (1.8*cos(x(1))*x(3)+(-0.8+cos(x(1))^2)*x(4))/(1+cos(x(1))^2)];
reset_lin = struct('A',[1,0;0,-0.75],'c',[0;0]);

% target
target = 1;

% transition
trans_lin = transition(guard_hyp,reset_lin,target);
trans_nonlin = transition(guard_ls,reset_nonlin,target);


% init location and display
loc_lin = location(inv_poly,trans_lin,sys_lin)

% init complex location and display
loc_nonlin = location(inv_ls,trans_nonlin,sys_nonlin)

% array of locations
[loc_lin;loc_nonlin]

% ------------------------------ END OF CODE ------------------------------
