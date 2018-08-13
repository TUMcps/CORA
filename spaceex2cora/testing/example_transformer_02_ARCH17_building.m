function example_transformer_02_ARCH17_building()
% example_linear_reach_ARCH17_building - example of linear reachability 
% analysis from the ARCH17 friendly competition (building example)
%
% Syntax:  
%    example_linear_reach_ARCH17_building
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
% Written:      09-February-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

load build

dim=length(A);

%set options --------------------------------------------------------------
R0 = interval([0.0002*ones(10,1); -0.0001; zeros(37,1)], [0.00025*ones(10,1); 0.0001; zeros(37,1)]);
options.x0=mid(R0); %initial state for simulation

%options.taylorTerms=4; %number of taylor terms for reachable sets
options.taylorTerms=6; %number of taylor terms for reachable sets
%options.zonotopeOrder=200; %zonotope order
options.zonotopeOrder=100; %zonotope order
options.originContained=0;
options.reductionTechnique='girard';
options.linAlg = 1;

uTrans=0.9;
options.u = uTrans;
options.uTrans=uTrans; %center of input set
options.U=zonotope([options.uTrans,0.1]); %input for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
buildingSys=linearSys('buildingSys',A,B); %initialize system
sha_build = SX2structHA('build_48.xml');
StructHA2file(sha_build,'buildingSX');
buildingSysSX = buildingSX();
%--------------------------------------------------------------------------

%compute reachable set using zonotopes
tic
%1st phase
options.R0=zonotope(R0); %initial state for reachability analysis
timeStep_1 = 0.002;
options.timeStep = timeStep_1;
%Rcont = reach(buildingSys, options);
%2ndphase
%options.R0 = Rcont{end};
% options.tStart=1; %start time
% options.tFinal=20; %final time
% timeStep_2 = 0.01;
% options.timeStep = timeStep_2;
% Rcont2 = reach(buildingSys, options);
% tComp = toc;
% disp(['computation time of reachable set: ',num2str(tComp)]);
% 
% %verification
% tic
% violation = 0;
% boundReached = 0;
% for i=1:length(Rcont)
%     x_25 = interval(project(Rcont{i},25));
%     if supremum(x_25) > 5.1e-3;
%         violation = 1;
%     end
%     if supremum(x_25) > 4e-3;
%         boundReached = 1;
%     end
% end
% for i=1:length(Rcont2)
%     x_25 = interval(project(Rcont2{i},25));
%     if supremum(x_25) > 5.1e-3;
%         violation = 1;
%     end
% end
% violation
% boundReached
% tVer = toc;
% disp(['computation time of verification: ',num2str(tVer)]);

options.tStart=0; %start time
options.tFinal=2; %final time


%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 10;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 20;
options.R0=zonotope(R0); %initial state for simulation results
%simRes = simulate_random(buildingSys, options, runs, fractionVertices, fractionInputVertices, inputChanges);
options.timeStep = 0.0001;

% simulate the origial system
[~,simRes.t,simRes.x] = simulate(buildingSys, options, 0, 2,options.x0,options);

% set further options that are necessary for hybrid automata
options.finalLoc = 2;
options.startLoc = 1;
options.Uloc = {zonotope(0)};
options.uLocTrans = {options.uTrans};

% simulate the system converted from spaceEX (= hybrid automata)
temp = simulate(buildingSysSX, options);
temp = get(temp,'result');
simResHybr = temp.simulation;

%plot results--------------------------------------------------------------
%plot results over time

figure;
hold on

% %plot time elapse
% for i=1:length(Rcont)
%     %get Uout 
%     Uout1 = interval(project(Rcont{i},25));
%     %obtain times
%     t1 = (i-1)*timeStep_1;
%     t2 = i*timeStep_1;
%     %generate plot areas as interval hulls
%     IH = interval([t1; infimum(Uout1)], [t2; supremum(Uout1)]);
% 
%     plotFilled(IH,[1 2],[.75 .75 .75],'EdgeColor','none');
% end
% 
% %plot time elapse, phase 2
% for i=1:length(Rcont2)
%     %get Uout 
%     Uout1 = interval(project(Rcont2{i},25));
%     %obtain times
%     t1 = options.tStart + (i-1)*timeStep_2;
%     t2 = options.tStart + i*timeStep_2;
%     %generate plot areas as interval hulls
%     IH = interval([t1; infimum(Uout1)], [t2; supremum(Uout1)]);
% 
%     plotFilled(IH,[1 2],[.75 .75 .75],'EdgeColor','none');
% end

% i = 1;
% iDim = 25;
% %plot time elapse
% while i<=length(Rcont)
%     %get Uout 
%     t1 = (i-1)*timeStep_1;
%     try
%         minVal = inf;
%         maxVal = -inf;
%         for k=1:10
%             Uout = interval(project(Rcont{i},iDim));
%             if infimum(Uout) < minVal
%                 minVal = infimum(Uout);
%             end
%             if supremum(Uout) > maxVal
%                 maxVal = supremum(Uout);
%             end
%             i = i + 1;
%         end
%     catch
%         minVal = infimum(interval(project(Rcont{i-1},iDim)));
%         maxVal = supremum(interval(project(Rcont{i-1},iDim)));
%     end
%     t2 = (i-1)*timeStep_1;
%     %generate plot areas as interval hulls
%     IH1 = interval([t1; minVal], [t2; maxVal]);
% 
%     plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
% end
% 
% i = 1;
% while i<=length(Rcont2)
%     %get Uout 
%     t1 = options.tStart + (i-1)*timeStep_2;
%     try
%         minVal = inf;
%         maxVal = -inf;
%         for k=1:5
%             Uout = interval(project(Rcont2{i},iDim));
%             if infimum(Uout) < minVal
%                 minVal = infimum(Uout);
%             end
%             if supremum(Uout) > maxVal
%                 maxVal = supremum(Uout);
%             end
%             i = i + 1;
%         end
%     catch
%         minVal = infimum(interval(project(Rcont2{i-1},iDim)));
%         maxVal = supremum(interval(project(Rcont2{i-1},iDim)));
%     end
%     t2 = options.tStart + (i-1)*timeStep_2;
%     %generate plot areas as interval hulls
%     IH1 = interval([t1; minVal], [t2; maxVal]);
% 
%     plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
% end

%plot simulation results
hold on
    plot(simRes.t,simRes.x(:,25),'Color',0*[1 1 1]);
    plot(simResHybr.t{1},simResHybr.x{1}(:,25),'r');

    temp = simRes.x(:,25)-simResHybr.x{1}(:,25);
    n = norm(temp);
    if(n > 1e-5)
        error('unit test failed, error=%i',n);
    else
        disp("n = " + string(n) + ", good job!");
    end
    
axis([0, 1, -6.5e-3, 6e-3])
%axis([0, 20, -6.5e-3, 6e-3])
xlabel('t');
ylabel('x_{25}');

% for plotRun=1:2
%     % plot different projections
%     if plotRun==1
%         projectedDimensions=[1 2];
%     elseif plotRun==2
%         projectedDimensions=[3 4];     
%     end 
%     
%     figure;
%     hold on
% 
%     %plot reachable sets 
%     for i=1:length(Rcont)
%         plotFilled(Rcont{i},projectedDimensions,[.8 .8 .8],'EdgeColor','none');
%     end
%     
%     %plot initial set
%     plot(options.R0,projectedDimensions,'w-','lineWidth',2);
%     
%     %plot simulation results      
%     for i=1:length(simRes.t)
%         plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
%     end
% 
%     %öabel plot
%     xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
%     ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
% end
%--------------------------------------------------------------------------

%------------- END OF CODE --------------
