function res = test_linearSysDT_observe_05_FRad_C()
% test_linearSysDT_observe_05_FRad_C - unit_test_function for guaranteed
%    state estimation of linear discrete-time systems (FRad-C method);
%    Checks the solution of the linearSysDT class for a four-dimensional 
%    vehicle model against an alternative implementation of the FRad-C
%    method published in [1].
%
% Syntax:
%    res = test_linearSysDT_observe_05_FRad_C
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Reference:
%    [1] C. Combastel. Zonotopes and Kalman observers:
%        Gain optimality under distinct uncertainty paradigms and
%        robust convergence. Automatica, 55:265-273, 2015.

% Authors:       Matthias Althoff
% Written:       25-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


%% Load side slip model
load slipEstimationModel_4D params options vehicle simRes

% set approach
options.alg = 'FRad-C';

%% observe
% standard implementation
estSet = observe(vehicle,params,options);
% alterantive implementation
estSet_alternative = setPropagationObserver_FRad_C_unitTest(vehicle,params,options);

%% compare whether enclosing hulls of the estimated sets match
% init
timeSteps = length(estSet.timePoint.set);
resPartial = zeros(timeSteps,1);
accuracy = 1e-8;
% loop over set
for iSet = 1:timeSteps
    % obtain interval hulls
    IH = interval(estSet.timePoint.set{iSet});
    IH_alternative = interval(estSet_alternative.timePoint.set{iSet});
    % check if slightly bloated versions enclose each other
    res_incl = isequal(IH,IH_alternative,1+accuracy);
    % check if simulation is enclosed
    res_sim = contains(IH,simRes.x{1}(iSet,:)');
    % combine results
    resPartial(iSet) = res_incl && res_sim;
end
% Are all sets matching?
res = all(resPartial);


% %% Plot results of estimated states
% % create time vector
% t = 0:options.timeStep:params.tFinal;
% % loop over all dimensions
% for iDim = 1:4
%     % create figure
%     figure
%     hold on
%     % plot estimated sets over time
%     h = plotOverTime(estSet,iDim,'FaceColor',[.6 .6 .6]);
%     % plot true state
%     p = plot(t(1:length(simRes.x{1}(:,iDim))),simRes.x{1}(:,iDim),'k.');
%     legend([h(end),p(end)],'Set of possible states','True state');
%     % axes labels
%     xlabel('t');
%     ylabel(['x_',num2str(iDim)]);
% end

        
% ------------------------------ END OF CODE ------------------------------
