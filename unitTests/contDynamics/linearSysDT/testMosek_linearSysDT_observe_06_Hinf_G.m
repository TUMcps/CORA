function res = testMosek_linearSysDT_observe_06_Hinf_G()
% testMosek_linearSysDT_observe_06_Hinf_G - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the 
% Hinf-G method and the Mosek 2021 solver.
%
% Checks the solution of the linearSysDT class for a four-dimensional 
% vehicle model against an alternative implementation of the Hinf-G method
% published in [1].
%
% Syntax:
%    res = testMosek_linearSysDT_observe_06_Hinf_G
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] W. Tang, Z. Wang, Y. Wang, T. Raissi, and Y. Shen.
%        Interval estimation methods for discrete-time linear time-
%        invariant systems. IEEE Transactions on Automatic Control,
%        64(11):4717-4724, 2019.

% Authors:       Matthias Althoff
% Written:       25-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


%% enable access to private function "observe_gain_Hinf"
path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSysDT','private','observe_gain_HinfG.m');
target = fullfile(path,'contDynamics','@linearSysDT','observe_gain_HinfG.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));


%% Load side slip model
load slipEstimationModel_4D params options vehicle simRes

% set approach
options.alg = 'Hinf-G';

%% observe
% standard implementation
estSet = observe(vehicle,params,options);
% alterantive implementation
options = params2options(params,options);
options.L = observe_gain_HinfG(vehicle,options);
estSet_alternative = intervalObserver_unitTest(vehicle,params,options);

%% compare whether enclosing hulls of the estimated sets match
% init
timeSteps = length(estSet.timeInterval.set);
resPartial = zeros(timeSteps,1);
accuracy = 1e-8;
% loop over set
for iSet = 1:timeSteps
    % obtain interval hulls
    IH = interval(estSet.timeInterval.set{iSet});
    IH_alternative = interval(estSet_alternative.timeInterval.set{iSet});
    % check if slightly bloated versions enclose each other
    res_cont = isequal(IH,IH_alternative,accuracy);
    % check if simulation is enclosed
    res_sim = contains(IH,simRes.x{1}(iSet,:)');
    % combine results
    resPartial(iSet) = res_cont && res_sim;
end
% Are all sets matching?
res = all(resPartial);

% revoke access to private function
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

% %% Plot results of estimated states
% % create time vector
% t = 0:options.timeStep:params.tFinal;
% % loop over all dimensions
% for iDim = 1:4
%     % create figure
%     figure
%     hold on
%     % plot estimated sets over time
%     h = plotOverTime(estSet,iDim,'FaceColor',[.6 .6 .6],'EdgeColor','none');
%     % plot true state
%     p = plot(t(1:length(simRes.x{1}(:,iDim))),simRes.x{1}(:,iDim),'k.');
%     legend([h(end),p(end)],'Set of possible states','True state');
%     % axes labels
%     xlabel('t');
%     ylabel(['x_',num2str(iDim)]);
% end

        
% ------------------------------ END OF CODE ------------------------------
