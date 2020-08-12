function res = test_reach_time()
% test_reach_time - unit_test_function of nonlinear reachability analysis
%    for following a reference trajectory
%
% Checks the solution of an autonomous car following a reference trajectory;
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
%    -

% Author:       Niklas Kochdumper
% Written:      23-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% Parameters --------------------------------------------------------------

dim_x = 8;
params.tFinal=4;
params.R0 = zonotope([[0; 0; 0; 22; 0 ; 0; -2.1854; 0],0.05*diag(ones(dim_x,1))]);
params.u = uTRansVec4CASreach();
params.u = params.u(:,1:400);
params.U = zonotope([0*params.u(:,1), 0.05*diag([ones(5,1);zeros(21,1)])]);


% Reachability settings ---------------------------------------------------

options.timeStep=0.01; %time step size for reachable set computation
options.taylorTerms=5; %number of taylor terms for reachable sets
options.zonotopeOrder=200; %zonotope order
options.maxError = ones(dim_x,1); % for comparison reasons
options.alg = 'lin';
options.tensorOrder = 2;
options.reductionInterval = Inf;


% System Dynamics ---------------------------------------------------------

vehicle = nonlinearSys(@vmodel_A_bicycle_linear_controlled); %initialize van-der-Pol oscillator


% Reachability Analysis --------------------------------------------------- 

Rset = reach(vehicle, params, options);


% Simulation --------------------------------------------------------------

simOpt.points = 100;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 1;
simOpt.inpChanges = 9;
simRes = simulateRandom(vehicle,params,simOpt);


% Numerical Evaluation ----------------------------------------------------

% check if end points are inside the final reachable set
R = Rset.timeInterval.set{end};
R = reduce(R,'girard',1);
R = halfspace(R);

for i = 1:length(simRes.x)
    temp = simRes.x{i}(end,:);
    res = all(R.halfspace.H*temp'<=R.halfspace.K);
    
    if ~res
       file_name = strcat('test_reach_time_', ...
                          datestr(now,'mm-dd-yyyy_HH-MM'));
                   
       file_path = fullfile(coraroot(), 'unitTests', 'failedTests', ...
                            file_name);
                  
       save(file_path, 'simRes')
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