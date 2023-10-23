function res = example_parallel_hybrid_02_lowPassFilter()
% example_parallel_hybrid_02_lowPassFilter - example for reachability of a
%    parallel hybrid automaton. The system consists of two piecewise linear
%    low-pass filters that are connected in series
%
% Syntax:
%    res = example_parallel_hybrid_02_lowPassFilter
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

% Authors:       Niklas Kochdumper
% Written:       06-July-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% System Dynamics ---------------------------------------------------------

pHA = lowpassFilter();


% Parameter ---------------------------------------------------------------

params.tFinal = 0.4;
params.startLoc = [1;3]; 
params.R0 = zonotope([0;0;0;0],diag([0.01,0.01,0.1,0.1]));  


% Reachability Settings ---------------------------------------------------

% settings for continuous reachability 
options.taylorTerms = 8; 
options.zonotopeOrder = 9; 
options.timeStep = 1e-04; 
 
% settings for hybrid systems
options.enclose = {'box'}; 
options.guardIntersect = 'conZonotope';
options.guardOrder = 3;


% Simulation --------------------------------------------------------------

simRes = simulateRandom(pHA,params);


% Reachability Analysis ---------------------------------------------------

R = reach(pHA,params,options);


% Visualization -----------------------------------------------------------

% plot filter k
projDim = {[1,2],[3,4]};
for k=1:2
    figure; box on; hold on;
    useCORAcolors("CORA:contDynamics")
    plot(R,projDim{k},'Unify',true,'UnifyTotalSets',5);
    updateColorIndex;
    plot(simRes,projDim{k});
    xlabel(['$x_{' num2str(projDim{k}(1)) '}$'],'interpreter','latex','FontSize',20);
    ylabel(['$x_{' num2str(projDim{k}(2)) '}$'],'interpreter','latex','FontSize',20);
    title(['Filter ' num2str(k)]);
end

res = true;

% ------------------------------ END OF CODE ------------------------------
