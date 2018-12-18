function res = test_nonlinear_reach_03_vanDerPol()
% test_nonlinear_reach_03_vanDerPol - unit_test_function of nonlinear reachability analysis
%
% Checks the solution of the nonlinearSys class for the van der Pol example;
% The settings are identical to the CDC'08 paper.
% It is checked whether the reachable set is enclosed in the initial set
% after a certain amount of time.
%
% Syntax:  
%    res = test_nonlinear_reach_03_vanDerPol()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      26-June-2009
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


%set options---------------------------------------------------------------
options.tStart=0; %start time
%options.tFinal=0.5; %final time
options.tFinal=6.76; %final time
options.x0=[1.4; 2.3]; %initial state for simulation

Z0{1}=zonotope([1.4 0.05 0; 2.3 0 0.05]); %initial state for reachability analysis
%Z0{1}=zonotope([1.4 0.15 0; 2.3 0 0.05]); %initial state for reachability analysis
%options.R0=zonotope([1.4 0.3 0; 2.3 0 0.05]); %initial state for reachability analysis
%options.R0=zonotope([1.4 0.6 0; 2.3 0 0.05]); %initial state for reachability analysis
options.R0 = zonotopeBundle(Z0);

options.timeStep=0.01; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=10; %zonotope order
options.intermediateOrder = 10;
options.errorOrder = 5;
options.polytopeOrder=1.5; %polytope order

options.advancedLinErrorComp = 1;
options.tensorOrder=2;
options.reductionTechnique='girard';
options.filterLength = [10, 5];
options.maxError=0.05*[1; 1];
%options.maxError=0.02*[1; 1];
options.reductionInterval=100;
%options.reductionInterval=50;

%uncertain inputs
options.uTrans = 0;
options.U = zonotope([0]); %input for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
vanderPol=nonlinearSys(2,1,@vanderPolEq,options); %initialize van-der-Pol oscillator
%--------------------------------------------------------------------------

% %simulate hybrid automaton
% simRes = simulate(vanderPol,options);

%plot(HA,'simulation',options);
      
tic
%compute reachable set
Rcont = reach(vanderPol, options);
toc

%obtain array of enclosing polytopes of last reachable set
for i=1:length(Rcont{end})
    if i==1
        Premain = enclosingPolytope(Rcont{end}{i},options);
    else
        Premain = Premain | enclosingPolytope(Rcont{end}{i},options);
    end
end

% % Visualization
% hold on
% for i = 1:20
%    for j = 1:length(Rcont{i})
%        plotFilled(Rcont{i}{j},[1,2],[.6 .6 .6],'EdgeColor','none');
%    end
% end
% plot(Premain,[1,2],'r');

%remove previous reachable sets
iStep = 1;
while ~isempty(Premain) && (iStep<=20)
    for iSet=1:length(Rcont{iStep})
        Preach = polytope(Rcont{iStep}{iSet});
        Premain = Premain\Preach;
    end
    iStep = iStep + 1;
end

%obtain result
if isempty(Premain)
    res = 1;
else
    res = 0;
end


% %saved result
% IH_saved = interval( ...
%            [2.616692934796753; 2.527784730195676; 2.371583574947125; 2.206019676085757; 2.063607715843610; 1.979956152875136], ...
%            [3.514750448595890; 3.317161121549810; 3.078098339981328; 2.833131900563080; 2.639101839639072; 2.564295199457605]);
%         
% %check if slightly bloated versions enclose each other
% res_1 = (IH <= enlarge(IH_saved,1+1e-8));
% res_2 = (IH_saved <= enlarge(IH,1+1e-8));
% 
% %final result
% res = res_1*res_2;


% %simulate
% stepsizeOptions = odeset('MaxStep',0.2*(options.tFinal-options.tStart));
% %generate overall options
% opt = odeset(stepsizeOptions);
% 
% %initialize
% runs=40;
% finalTime = options.tFinal;
% 
% for i=1:runs
% 
%     %set initial state, input
%     if i<=30
%         options.x0=randPointExtreme(Z0{1}); %initial state for simulation
%     else
%         options.x0=randPoint(Z0{1}); %initial state for simulation
%     end
% 
%     %set input
%     if i<=8
%         options.u=randPointExtreme(options.U)+options.uTrans; %input for simulation
%     else
%         options.u=randPoint(options.U)+options.uTrans; %input for simulation
%     end
% 
%     %simulate hybrid automaton
%     [vanderPol,t{i},x{i}] = simulate(vanderPol,options,options.tStart,options.tFinal,options.x0,opt); 
% end
% 
% 
% plotOrder = 20;
% 
% %plot dynamic variables
% for plotRun=1:1
% 
%     
%     if plotRun==1
%         projectedDimensions=[1 2];
%     end 
% 
%     figure;
%     hold on
% 
%     %plot reachable sets of zonotope; mode 1
%     for i=1:length(Rcont)
%         for j=1:length(Rcont{i})
%             Zproj = reduce(Rcont{i}{j},'girard',plotOrder);
%             %plot(Zproj,[1 2],'lightgray');
%             plot(Zproj,[1 2],'b');
%         end
%     end
%     
%   
%     %plot simulation results      
%     for i=1:length(x)
%         plot(x{i}(:,1),x{i}(:,2),'k');
%     end
%     
%     %plot initial set
%     %plot(options.R0,projectedDimensions,'whiteEdge');
% 
%     xlabel('x_1');
%     ylabel('x_2');
% end


%------------- END OF CODE --------------