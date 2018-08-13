function completed = example_linear_reach_05_discreteTime()
% example_linear_reach_05_discreteTime - example for linear reachability 
% analysis using discrete time
%
% Syntax:  
%    completed = example_linear_reach_05_discreteTime()
%
% Inputs:
%    no
%
% Outputs:
%    completed - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      31-July-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=5; %final time
options.R0=zonotope([ones(5,1),0.1*eye(5)]); %initial state for reachability analysis
options.timeStep=0.04; %time step size for reachable set computation
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
A=[0.8, -0.4, 0, 0, 0; ...
   0.4, 0.8, 0, 0, 0; ...
   0, 0, 0.7, 0.07, 0; ...
   0, 0, -0.07, 0.7, 0; ...
   0, 0, 0, 0, 0.8]; % system matrix
U = zonotope([zeros(5,1), 0.1*diag([0.2, 0.5, 0.2, 0.5, 0.5])]); % uncertain input
%--------------------------------------------------------------------------

%compute reachable set using zonotopes
tic
t = options.tStart;
i = 1;
R{1} = options.R0;
while t < options.tFinal %final time not yet reached
    R_next = A*R{i} + U;
    R{i+1} = reduce(R_next,'girard', 10);
    i = i+1;
    t = t+options.timeStep;
end
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%plot results--------------------------------------------------------------
for plotRun=1:2
    % plot different projections
    if plotRun==1
        projectedDimensions=[1 2];
    elseif plotRun==2
        projectedDimensions=[3 4];     
    end 
    
    figure;
    hold on

    %plot reachable sets 
    for i=1:length(R)
        plotFilled(R{i},projectedDimensions,[.8 .8 .8],'EdgeColor','none');
    end
    
    %plot initial set
    plot(options.R0,projectedDimensions,'w-','lineWidth',2);

    %label plot
    xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
    ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
end
%--------------------------------------------------------------------------

%example completed
completed = 1;


%------------- END OF CODE --------------
