function res = test_reach_time()
% test_reach_time - unit_test_function of nonlinear reachability analysis
% for following a reference trajectory
%
% Checks the solution of an autonomous car following a reference
% trajectory;
% It is checked whether the final reachable set encloses the end points of
% the simulated trajectories
%
% Syntax:  
%    res = test_reach_time()
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
% Author:       Niklas Kochdumper
% Written:      23-March-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


dim = 8;

options.maxError = ones(dim,1); % for comparison reasons


%set options --------------------------------------------------------------
options.tStart=0; %start time
%options.tFinal=3.99; %final time
options.tFinal=0.1; %final time
options.x0=[0; 0; 0; 22; 0 ; 0; -2.1854; 0]; %initial state for simulation
options.R0 = zonotope([options.x0, 0.05*diag([1, 1, 1, 1, 1, 1, 1, 1])]); %initial state for reachability analysiszonotope([options.x0, diag([0.20, 0.20])]); %max for 3rd order
options.timeStep=0.01; %time step size for reachable set computation
options.taylorTerms=5; %number of taylor terms for reachable sets
options.zonotopeOrder=200; %zonotope order

options.uTransVec = uTRansVec4CASreach();
options.u = 0;
options.U = zonotope([0*options.uTransVec(:,1), 0.05*diag([ones(5,1);zeros(21,1)])]);

options.advancedLinErrorComp = 0;
options.tensorOrder = 1;
options.reductionInterval = inf;
options.reductionTechnique = 'girard';
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
vehicle = nonlinearSys(8,26,@vmodel_A_bicycle_linear_controlled,options); %initialize van-der-Pol oscillator
%--------------------------------------------------------------------------

%compute reachable set 
[~,Rset] = reach(vehicle, options);

% simulate the system
options.uTransVec = options.uTransVec(:,1:10);

simRes = simulate_random(vehicle,options,100,0.5,1,9);

% check if end points are inside the final reachable set
R = Rset{end}{1}.set;
R = reduce(R,'girard',1);
R = halfspace(R);

for i = 1:length(simRes.x)
    temp = simRes.x{i}(end,:);
    res = all(R.halfspace.H*temp'<=R.halfspace.K);
    if res == 0
       break; 
    end
end

% plot the result
% counter = 1;
% for j = 1:4
%     figure
%     plot(Rset{end}{1}.set,[counter,counter+1],'r');
%     hold on
%     for i = 1:length(simRes.x)
%        temp = simRes.x{i}(end,:);
%        plot(temp(counter),temp(counter+1),'.k');
%     end
%     counter = counter + 2;
% end