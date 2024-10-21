function res = example_nonlinear_reachInner_01_brusselator
% example_nonlinear_reachInner_01_brusselator - example for the computation
%    of an inner approximation of the reachable set for nonlinear dynamics,
%    using the Minkowski difference approach from [1]
%
% Syntax:
%    res = example_nonlinear_reachInner_01_brusselator
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References:
%    [1] M. Wetzlinger, A. Kulmburg, and M. Althoff. "Inner approximations
%        of reachable sets for nonlinear systems using the Minkowski
%        difference". IEEE Control Systems Letters, 2024.

% Authors:       Mark Wetzlinger
% Written:       17-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------

brusselator = @(x,u) [1-2*x(1) + 3/2 * x(1)^2*x(2); ...
                      x(1)-3/2*x(1)^2*x(2)];
sys = nonlinearSys('brusselator',brusselator);


% Parameters --------------------------------------------------------------

params.tFinal = 3;
R0 = interval([0.9;0],[1;0.1]);


% Reachability Settings ---------------------------------------------------

options_inner_Minkdiff.algInner = 'minkdiff';
options_inner_Minkdiff.timeStep = 0.01;
options_inner_Minkdiff.compOutputSet = false;
options_inner_Minkdiff.tensorOrder = 2;


% Reachability Analysis ---------------------------------------------------

params.R0 = polytope(R0);
Rin = reachInner(sys,params,options_inner_Minkdiff);


% Visualization -----------------------------------------------------------

figure; hold on; box on;
xlabel('x_1'); ylabel('x_2');

% inner approximation
useCORAcolors("CORA:contDynamics");
plot(Rin,[1,2],'DisplayName','Inner approximation');
plot(Rin.R0,[1,2],'DisplayName','Initial set');

% completed
res = true;

% ------------------------------ END OF CODE ------------------------------
