function res = testLongDuration_linearSysDT_observe_04_FRad_A()
% testLongDuration_linearSysDT_observe_04_FRad_A - unit_test_function for guaranteed
% state estimation of linear discrete-time systems using the FRad-A method.
%
% Checks the solution of the linearSysDT class for a four-dimensional 
% vehicle model against an alternative implementation of the FRad-A method
% published in [1].
%
% Syntax:  
%    res = testLongDuration_linearSysDT_observe_04_FRad_A
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed
%        state estimation by zonotopes. Automatica, 41(6):1035-1043,
%        2005.
%
% Example: 
%    -
 
% Author:       Matthias Althoff
% Written:      25-Feb-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%% Load side slip model
load slipEstimationModel_4D params options vehicle simRes

% set approach
options.alg = 'FRad-A';

%% observe
% standard implementation
estSet = observe(vehicle,params,options);
% alterantive implementation
estSet_alternative = stripBasedObserver_unitTest(vehicle,params,options);

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
    res_encl = (IH <= enlarge(IH_alternative,1+accuracy));
    res_incl = (IH_alternative <= enlarge(IH,1+accuracy));
    % check if simulation is enclosed
    res_sim = in(IH,simRes.x{1}(iSet,:)');
    % combine results
    resPartial(iSet) = res_encl*res_incl*res_sim;
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
%     h = plotOverTime(estSet,iDim,'FaceColor',[.6 .6 .6],'EdgeColor','none');
%     % plot true state
%     p = plot(t(1:length(simRes.x{1}(:,iDim))),simRes.x{1}(:,iDim),'k.');
%     legend([h(end),p(end)],'Set of possible states','True state');
%     % axes labels
%     xlabel('t');
%     ylabel(['x_',num2str(iDim)]);
% end

        
%------------- END OF CODE --------------
