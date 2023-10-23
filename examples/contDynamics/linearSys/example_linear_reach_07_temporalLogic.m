function res = example_linear_reach_07_temporalLogic
% example_linear_reach_07_temporalLogic - example for model checking 
%    reachable set against temporal logic specifications according to [1]
%
% Syntax:
%    res = example_linear_reach_07_temporalLogic
%
% Inputs:
%    ---
%
% Outputs:
%    res - true/false
%
% References:
%    [1] H. Roehm et al. "STL Model Checking of Continuous and Hybrid
%        Systems", International Symposium on Automated Technology for 
%        Verification and Analysis, pp. 412-427, 2016.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       16-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.R0 = zonotope([0;-1],diag([0.1,0.1]));
params.tFinal = 2;


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.25;
options.zonotopeOrder = 10;
options.taylorTerms = 4;


% System Dynamics ---------------------------------------------------------

A = [0 -1; 1 0];
B = [0;0];

sys = linearSys(A,B);


% Specification -----------------------------------------------------------

x = stl('x',2);
eq = until(x(2) < -0.5,x(1) > 0.5,interval(0,2));

spec = specification(eq,'logic');


% Reachability Analysis ---------------------------------------------------

tic
[R,res] = reach(sys, params, options, spec);
tComp = toc;

disp(['computation time of reachable set: ',num2str(tComp)]);
if res
    disp('specification satisfied!');
else
    disp('specification violated!');
end


% Visualization -----------------------------------------------------------

% formatting
figure; hold on; box on;
xlim([-1,2]); ylim([-2,1]);
xlabel('x_1'); ylabel('x_2');

% plot halfspaces for the specifications
hs1 = convert2set(eq.lhs);
hs2 = convert2set(eq.rhs);

plot(hs1,[1,2],'FaceColor',colorblind('r'),'FaceAlpha',0.5);
plot(hs2,[1,2],'FaceColor',colorblind('b'),'FaceAlpha',0.5);

% plot reachable set
useCORAcolors("CORA:contDynamics")
plot(R,[1,2],'DisplayName','Reachable set');

% plot initial set
plot(R.R0,[1,2],'DisplayName','Initial set');

% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
