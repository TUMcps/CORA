function res = example_nonlinear_reachInner_02_biologicalModel
% example_nonlinear_reachInner_02_biologicalModel - example for the
%    computation of an inner approximation of the reachable set for
%    nonlinear dynamics, using the approaches in [1] and [2]
%
% Syntax:
%    res = example_nonlinear_reachInner_02_biologicalModel
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

sys = nonlinearSys('biologicalModel',@biologicalModel);


% Parameters --------------------------------------------------------------

params.tFinal = 0.2;
R0 = interval(0.99*ones(7,1),1.01*ones(7,1));


% Reachability Settings ---------------------------------------------------

options_outer.alg = 'lin-adaptive';

% 1. Minkdiff algorithm
options_inner_Minkdiff.algInner = 'minkdiff';
options_inner_Minkdiff.timeStep = 0.01;
options_inner_Minkdiff.tensorOrder = 2;
options_inner_Minkdiff.compOutputSet = false;

% 2. polyZonotope/scale (from CodeOcean capsule)
options_inner_scale.algInner = 'scale';
options_inner_scale.timeStep = 0.0001;
options_inner_scale.taylorTerms = 10;
options_inner_scale.zonotopeOrder = 30;
options_inner_scale.intermediateOrder = 10;
options_inner_scale.errorOrder = 5;
options_inner_scale.splits = 8;
options_inner_scale.iter = 2;
options_inner_scale.orderInner = 2;
options_inner_scale.scaleFac = 0.9;


% Reachability Analysis ---------------------------------------------------

params.R0 = zonotope(R0);
Rout = reach(sys,params,options_outer);

params.R0 = polytope(R0);
Rin_Minkdiff = reachInner(sys,params,options_inner_Minkdiff);

params.R0 = polyZonotope(R0);
Rin_scale = reachInner(sys,params,options_inner_scale);


% Visualization -----------------------------------------------------------

figure;
projDims = {[1,2],[3,4],[5,6],[1,7],[2,3],[4,5],[6,7],[1,3],[2,4]};

for p=1:length(projDims)

    % labels
    subplot(3,3,p); hold on; box on;
    xlabel("x_" + projDims{p}(1)); ylabel("x_" + projDims{p}(2));
    useCORAcolors('CORA:contDynamics',3);
    
    % outer approximation 
    plot(Rout,projDims{p},'DisplayName','Outer approximation');
    
    % inner approximation
    plot(Rin_Minkdiff,projDims{p},'DisplayName','Inner approximation (Minkowski difference)');
    plot(Rin_scale,projDims{p},'DisplayName','Inner approximation (Scaling)');

    % initial set
    plot(Rout.R0,projDims{p},'DisplayName','Initial set');

end

% completed
res = true;

% ------------------------------ END OF CODE ------------------------------
