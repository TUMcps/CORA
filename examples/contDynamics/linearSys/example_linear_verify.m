function res = example_linear_verify
% example_linear_verify - test function for the novel verification
%    algorithm based on rigorous Hausdorff distance to exact solution
%
% Syntax:  
%    res = example_linear_verify
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      16-August-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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
hs = halfspace([0 -1],-12.1);
spec = specification(hs,'unsafeSet');

P = mptPolytope([-1 0;0 -1;1 2],[0;0;4]) + [-3.8;5.833];
spec = add(spec,specification(P,'unsafeSet'));

% safe sets
P = mptPolytope([1 0;-1 0;0 1;0 -1;1 -1;-1 1], ...
                        [12; 7; 12.2; -2.05; 1.2; 13]);
spec = add(spec,specification(P,'safeSet'));


% Verification ------------------------------------------------------------

tic;
[res,R] = verify(sys,params,spec);
tComp = toc;
disp(['Computation time: ',num2str(tComp),'s']);
disp(['Specifications satisfied? ' num2str(res)]);


% Visualization -----------------------------------------------------------

% specifications
figure; hold on; box on;
% goal set
h_safe = plot(spec(3).set,[1,2],'FaceColor',colorblind('gray'),'EdgeColor','k');
% unsafe sets
h_unsafe = plot(spec(1).set,[1,2],'FaceColor',colorblind('r'));
plot(spec(2).set,[1,2],'FaceColor',colorblind('r'));

% outer-approximation reachable set
h_outer = plot(R,[1,2]);

% initial set
h_init = plot(params.R0,[1,2],'FaceColor','w','EdgeColor','k');

% legend
legend([h_safe,h_unsafe,h_outer,h_init],'Safe set','Unsafe sets',...
    'Outer-approximation','Initial set','location','southeast');

%------------- END OF CODE --------------