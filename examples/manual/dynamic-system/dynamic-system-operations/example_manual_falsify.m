function example_manual_falsify()
% example_manual_falsify - example from the manual demonstrating the 
% falsify operation as defined in the manual
%
% Syntax:
%   example_manual_falsify()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Niklas Kochdumper
% Written:       15-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system dynamics
sys = linearSys([-0.7 -2;2 -0.7],[1;1],[-2;-1]);

% parameter
params.tFinal = 5;
params.R0 = zonotope(interval([2;2],[2.5;2.5]));
params.U = zonotope(interval(-0.1,0.1));

% specification
I = interval([0.7;-1],[1.4;0]);
spec = specification(I,'unsafeSet');

% falsification
[res,fals] = falsify(sys,params,spec);

% plot --------------------------------------------------------------------

figure; hold on;
plot(params.R0,[1,2],'FaceColor','w','EdgeColor','k');
plot(spec.set);
useCORAcolors("CORA:contDynamics")
plot(fals.traj);

xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
