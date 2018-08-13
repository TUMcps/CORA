function example_linearParam_reach_ARCH17_platoon_unbounded()
% example_linearParam_reach_ARCH17_platoon -  example of linear reachability 
% analysis from the ARCH17 friendly competition (platoon example); the
% linear dynamics can switch arbitrarily
%
% Syntax:  
%    example_linearParam_reach_ARCH17_platoon
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

% controlled system
A_c = [...
        0    1.0000         0         0         0         0         0         0         0;...
        0         0   -1.0000         0         0         0         0         0         0;...
        1.6050    4.8680   -3.5754   -0.8198    0.4270   -0.0450   -0.1942    0.3626   -0.0946;...
        0         0         0         0    1.0000         0         0         0         0;...
        0         0    1.0000         0         0   -1.0000         0         0         0;...
        0.8718    3.8140   -0.0754    1.1936    3.6258   -3.2396   -0.5950    0.1294   -0.0796;...
        0         0         0         0         0         0         0    1.0000         0;...
        0         0         0         0         0    1.0000         0         0   -1.0000;...
        0.7132    3.5730   -0.0964    0.8472    3.2568   -0.0876    1.2726    3.0720   -3.1356 ]; 

% uncontrolled system
A_n = [...
        0    1.0000         0         0         0         0         0         0         0;...
        0         0   -1.0000         0         0         0         0         0         0;...
        1.6050    4.8680   -3.5754         0         0         0         0         0         0;...
        0         0         0         0    1.0000         0         0         0         0;...
        0         0    1.0000         0         0   -1.0000         0         0         0;...
        0         0         0    1.1936    3.6258   -3.2396         0         0         0;...
        0         0         0         0         0         0         0    1.0000         0;...
        0         0         0         0         0    1.0000         0         0   -1.0000;...
        0.7132    3.5730   -0.0964    0.8472    3.2568   -0.0876    1.2726    3.0720   -3.1356 ];   
    
    
B = [0 ; 1; 0; 0; 0; 0; 0; 0; 0 ];
        
%build zonotpe matrix
A_mid = 0.5*(A_c + A_n);
A_gen{1} = 0.5*(A_c - A_n);
matZ_A = matZonotope(A_mid, A_gen);

%get dimension
dim=matZ_A.dim;

%initial set
R0 = zonotope(zeros(dim,1));

%initial set
options.x0=center(R0); %initial state for simulation
options.R0=R0; %initial state for reachability analysis

%inputs
u = interval(-9,1); %[-9,1]
U = B*zonotope(u);
options.uTrans=0;
options.U=U+(-options.uTrans); %input for reachability analysis

%other
options.tStart=0; %start time
options.tFinal=50; %final time; prev 80
%options.intermediateOrder = 3;
options.intermediateOrder = -1; % did not really help
options.reductionTechnique='girard';

options.originContained = 1;
options.timeStep = 0.01;
%options.timeStep = 0.005; % did not really help
options.eAt = expm(matZ_A.center*options.timeStep);

%options.zonotopeOrder=200; %zonotope order
options.zonotopeOrder=400; %zonotope order
%options.zonotopeOrder=1000; %zonotope order; helped a little
options.taylorTerms=3;
%options.taylorTerms=2;

%time step
r = options.timeStep;
maxOrder=options.taylorTerms;

%instantiate linear dynamics with constant parameters
linSys  = linVarSys(matZ_A, eye(dim), r, maxOrder);

%reachable set computations
tic
[Rcont] = reach(linSys, options);

%invariant computation
Rlast = reduce(Rcont{end},'girard', 400);
options.R0 = enlarge(Rlast,1.01);
options.zonotopeOrder=800;
[RcontInv,tadd] = reachInv(linSys, options, options.R0);

tComp = toc;
disp(['computation time for platoon: ',num2str(tComp)]);

%verification
tic
violation30 = 0;
violation42 = 0;
violation50 = 0;
for i=1:length(Rcont)
    x_proj = interval(project(Rcont{i},[1,4,7]));
    if any(infimum(x_proj) < -30)
        violation30 = 1;
    end
    if any(infimum(x_proj) < -42)
        violation42 = 1;
    end
    if any(infimum(x_proj) < -50)
        violation50 = 1;
    end
end
for i=1:length(RcontInv)
    x_proj = interval(project(Rcont{i},[1,4,7]));
    if any(infimum(x_proj) < -30)
        violation30 = 1;
    end
    if any(infimum(x_proj) < -42)
        violation42 = 1;
    end
    if any(infimum(x_proj) < -50)
        violation50 = 1;
    end
end
violation30 
violation42
violation50 
tVer = toc;
disp(['computation time of verification: ',num2str(tVer)]);


%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
options.R0 = R0;
options.tFinal=options.tFinal + tadd; %final time
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(linSys, options, runs, fractionVertices, fractionInputVertices, inputChanges);

%plot results--------------------------------------------------------------
for plotRun=1:0
    % plot different projections
    if plotRun==1
        projectedDimensions=[1 2];
    elseif plotRun==2
        projectedDimensions=[3 4];
    elseif plotRun==3
        projectedDimensions=[5 6];
    elseif plotRun==4
        projectedDimensions=[7 8];
    else
        projectedDimensions=[8 9];
    end
    
    figure;
    hold on

    for i=1:length(Rcont)
        Zproj = project(Rcont{i},projectedDimensions);
        Zproj = reduce(Zproj,'girard',3);
        plotFilled(Zproj,[1 2],[.75 .75 .75],'EdgeColor','none');
    end

    %plot initial set
    plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');
       
    %plot simulation results      
    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
    end

    %öabel plot
    xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
    ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
end

%plot results over time
for iDim = 1:1

    figure;
    hold on

    %plot time elapse
    for i=1:length(Rcont)
        %get Uout 
        Uout1 = interval(project(Rcont{i},iDim));
        %obtain times
        t1 = (i-1)*options.timeStep;
        t2 = i*options.timeStep;
        %generate plot areas as interval hulls
        IH1 = interval([t1; infimum(Uout1)], [t2; supremum(Uout1)]);

        plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
    end
    for i=1:length(RcontInv)
        %get Uout 
        Uout1 = interval(project(RcontInv{i},iDim));
        %obtain times
        t1 = options.tFinal - tadd + (i-1)*options.timeStep;
        t2 = options.tFinal - tadd + i*options.timeStep;
        %generate plot areas as interval hulls
        IH1 = interval([t1; infimum(Uout1)], [t2; supremum(Uout1)]);

        plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
    end

    %plot simulation results
    for i=1:(length(simRes.t))
        plot(simRes.t{i},simRes.x{i}(:,iDim),'Color',0*[1 1 1]);
    end

end

%--------------------------------------------------------------------------


function [Rcont,t] = reachInv(obj,options,inv)


%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %time step
    r = options.timeStep;
    %compute initial state factor
    options.factor(i)= r^(i)/factorial(i);    
end
%possibility for updating time step
%options.timeStep = options.timeFactor*norm(A^2,inf)^(-0.5); 


%if a trajectory should be tracked
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
end

%initialize reachable set computations
[Rnext, options] = initReach(obj, options.R0, options);

%while final time is not reached
t=options.tStart;
iSet=1;
notInInv=1;

while notInInv
    
    %save reachable set in cell structure
    Rcont{iSet} = Rnext.ti; 
    Rcont_tp{iSet} = Rnext.tp; 
    
    %increment time and set counter
    t = t+options.timeStep;
    iSet = iSet+1; 
    options.t=t;
    if isfield(options,'verbose') && options.verbose 
        disp(t); %plot time
    end
    
    %if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,iSet);
    end
    
    %compute next reachable set
    [Rnext,options]=post(obj,Rnext,options);
    
    if inViaProj(Rnext.ti, inv)
        notInInv=0;
    end
    
end



%------------- END OF CODE --------------