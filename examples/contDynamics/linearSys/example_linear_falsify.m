function res = example_linear_falsify
% example_linear_falsify - example for safety falsification
%
% Syntax:
%    res = example_linear_falsify
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
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       15-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------

A = [-0.7 -2; 2 -0.7];
B = 1;
sys = linearSys('sys',A,B);


% Parameters --------------------------------------------------------------

% time horizon
params.tFinal = 1;

% initial set
params.R0 = zonotope(10*ones(2,1),0.5*diag(ones(2,1)));

% uncertain inputs
params.U = zonotope([zeros(2,1),0.05*eye(2)]);

params.u = repmat([0.1,-0.2,0.3,-0.4],2,1);
params.tu = [0,0.5,0.7,0.9];


% Specifications ----------------------------------------------------------

% unsafe sets
hs = polytope([0 -1],-12.1);
spec = specification(hs,'unsafeSet');

P = polytope([-1 0;0 -1;1 2],[0;0;4]) + [-3.8;5.833];
spec = add(spec,specification(P,'unsafeSet'));

% safe sets
P = polytope([1 0;-1 0;0 1;0 -1;1 -1;-1 1], ...
                        [12; 6.8; 12.2; -2.05; 1.2; 13]);
spec = add(spec,specification(P,'safeSet'));


% Falsification -----------------------------------------------------------

timerVal = tic;
[res,fals] = falsify(sys,params,spec);
tComp = toc(timerVal);
disp(['Computation time: ',num2str(tComp),'s']);
disp(['Falsification successful? ' num2str(res)]);


% Visualization -----------------------------------------------------------

% specifications
figure; hold on; box on;
plot(spec(3),[1,2], 'DisplayName','Safe set');
plot(spec(1:2),[1,2], 'DisplayName', 'Unsafe sets');

% initial set
plot(params.R0,[1,2],'FaceColor','w','EdgeColor','k', ...
                                        'DisplayName','Initial set');

% falsifying trajectory
useCORAcolors("CORA:contDynamics")
plot(fals.traj,[1,2],'DisplayName','Falsifying trajectory');


% legend
legend('location','southeast');

% ------------------------------ END OF CODE ------------------------------
