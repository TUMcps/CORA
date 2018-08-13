function completed = example_reachNonLinearDrumSystem()
% As realized in the revised paper of the IEEE PST elguindy-2017-PST
% Formal analyis of Drum-boiler units to maximize the load-following
% capabliiltes of power plants
%
% Syntax:  
%    ReachNonLinearDrumSystem
%
% Inputs:
%    no
%
% Outputs:
%    no 
%
% Example: 
% 
% 
% Author:        Ahmed El-Guindy
% Written:       11-November-2015
% Last update:   21-April-2018 (substantially implified, MA)
% Last revision: 14-July-2017

%------------- BEGIN CODE --------------

% Number of state variables of the drum
dim = 12;

% Set options -------------------------------------------------------------
options.tStart   = 0;   % Start time
options.timeStep = 1;   % Time step size for reachable set computation
options.tFinal   = 300; % Final time

% Time-horizon (t = 5 mins for the TSOs as per secondary frequency control 
% in Germany
TSim = (1:1:options.tFinal)';

% Initial set for a change of 40MW from 110MW towards 70MW
options.x0=[5.5e5 ;19.06; 0.01027 ;3.687; 967.7; 9.241e5;22e6;39.88;39.88;14.44;6.382;6.382]; 

%Initial set for a change of 40MW from 70MW towards 110MW
%options.x0=[5.5e5 ;18.06; 0.01427 ;2.687; 867.7; 8.241e5;18e6;39.88;39.88;20.44;6.382;6.382];

% Mapping of the reachable set to obtain the water level of the drum
C = [    1e-5       0          0           0         0 0 0       0 0 0 0 0;
         0          0          0           0         0 0 4.75e-6 0 0 0 0 0;
         0          0          0           0         0 0 0       0 0 0 1 0;
        -4.7523e-04 6.8027e+01 1.9778e+04 6.8027e+01 0 0 0       0 0 0 0 0 ]; 


% Uncertain initial sets
Bound(1,1) = 0.1e5;  Bound(2,2) = 0.5;Bound(3,3) = 0.001;Bound(4,4) = 0.5;
Bound(7,7) = 0.00001;Bound(12,12)=0;  Bound(7,7)=0.1e6;  Bound(8,8)=2;
Bound(10,10)=1;      Bound(11,11)=1;

% Initial sets for reachability analysis
options.R0 = zonotope([options.x0,Bound]); 
R0Level    = C*options.R0 + [ 0;0;0 ;(-1.4778e3)];

options.zonotopeOrder = 500; % Zonotope order
options.polytopeOrder = 10;  % Polytope order
options.taylorTerms   = 20;


% Initial set for the input variables
options.uTrans = [70; 0; 0];

% Uncertainty
Bound_u(1,1)  = 0; 
Bound_u(2,2)  = 0.01;
Bound_u(3,3)  = 5; 

options.U    = zonotope([zeros(3,1),Bound_u]);
options.originContained = 0;

% Linearization options
options.tensorOrder          = 2;
options.advancedLinErrorComp = 1;
options.intermediateOrder    = 3;
options.reductionTechnique   = 'girard';

options.errorOrder           = 1.5;
options.oldError             = zeros(dim,1);
options.maxError             = 1.0e+02 * ones(dim,1);
options.reductionInterval    = inf;

%--------------------------------------------------------------------------
%specify continuous dynamics-----------------------------------------------
DrumSys = nonlinearSys(12,3,@DrumModel,options);
save DrumSys DrumSys
load DrumSys
%--------------------------------------------------------------------------

% Compute reachable set 
tic
Rcont = reach(DrumSys, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


RcontY{length(Rcont)} = 0;
% Mapping to get the rechable set of the water level 
for i = 1:length(Rcont)
    RcontY{i} = C*Rcont{i}{1} + [ 0;0; 0;(-1.4778e3)];
end

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(DrumSys, options, runs, fractionVertices, fractionInputVertices, inputChanges);
% obtain level from states
for iRun = 1:length(simRes.x)
    level{iRun} = -4.7523e-04*simRes.x{iRun}(:,1) + 6.8027e+01*simRes.x{iRun}(:,2) + 6.8027e+01*simRes.x{iRun}(:,4) + 1.9778e+04*simRes.x{iRun}(:,3) - 1.4778e3;
end 


%% Plotting
%--------------------------------------------------------------------------
figure
subplot(2,2,1)
hold on;  


% Power demand vs level
for i=1:length(RcontY)    
    ZprojY = project(RcontY{i},[2 4]);
    ZprojY = reduce(ZprojY,'girard',10);
    plotFilled(ZprojY,[1 2],[.8 .8 .8],'EdgeColor','none');  
end
plotFilled(R0Level,[2 4],'w','EdgeColor','k');
for i=1:length(simRes.t)
    plot(simRes.x{i}(:,7)*4.75e-6,level{i},'Color',0*[1 1 1]);
end

subplot(2,2,2)
hold on; 

% Pressure vs level
for i=1:length(RcontY)    
    ZprojY = project(RcontY{i},[1 4]);
    ZprojY = reduce(ZprojY,'girard',10);
    plotFilled(ZprojY,[1 2],[.8 .8 .8],'EdgeColor','none');   
end
plotFilled(R0Level,[1 4],'w','EdgeColor','k');

for i=1:length(simRes.t)
    plot(simRes.x{i}(:,1)*1e-5,level{i},'Color',0*[1 1 1]);
end

subplot(2,2,3)
hold on; 

% Power demand vs Pressure
for i=1:length(RcontY)    
    ZprojY = project(RcontY{i},[2 1]);
    ZprojY = reduce(ZprojY,'girard',10);
    plotFilled(ZprojY,[1 2],[.8 .8 .8],'EdgeColor','none');   
end
plotFilled(R0Level,[2 1],'w','EdgeColor','k');
for i=1:length(simRes.t)
    plot(simRes.x{i}(:,7)*4.75e-6,simRes.x{i}(:,1)*1e-5,'Color',0*[1 1 1]);
end

projectedDimensions=[3 4];
subplot(2,2,4)
hold on;

% Steam quality vs Steam volume
    for i=1:length(Rcont)
        Zproj = project(Rcont{i}{1},projectedDimensions);
        Zproj = reduce(Zproj,'girard',10);
        plotFilled(Zproj,[1 2],[.8 .8 .8],'EdgeColor','none');        
    end

    plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');

    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
    end
    
    xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
    ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
    
    
%example completed
completed = 1;   

%------------- END OF CODE --------------