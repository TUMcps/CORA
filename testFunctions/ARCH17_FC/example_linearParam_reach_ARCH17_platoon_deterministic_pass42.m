function example_linearParam_reach_ARCH17_platoon_deterministic_pass42()
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
% Written:      05-April-2017
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
        
%get dimension
dim=length(A_c);

%initial set
Rinit = zonotope(zeros(dim,1));
Rtp{1} = Rinit;

%initial set
options.x0=center(Rtp{1}); %initial state for simulation

%inputs
u = interval(-9,1); %[-9,1]
U = B*zonotope(u);
options.uTrans=center(U);
options.U=U+(-options.uTrans); %input for reachability analysis

%other
options.tStart=0; %start time
options.tFinal=5; %final time
options.reductionTechnique='girard';
options.originContained = 0;
options.timeStep = 0.02; 

options.zonotopeOrder=20; %zonotope order
options.polytopeOrder=3; %polytope order
options.taylorTerms=4;
options.linAlg = 0;

%instantiate linear dynamics with constant parameters
linSys_c  = linearSys('c',A_c, eye(dim));
linSys_n  = linearSys('n',A_n, eye(dim));

%reachable set computations
tic
for i=1:2
    options.R0 = Rtp{end};
    [Rcont{2*i-1},Rtp] = reach(linSys_c, options);
    options.R0 = Rtp{end};
    [Rcont{2*i},Rtp] = reach(linSys_n, options);
end
tComp = toc;
disp(['computation time for platoon: ',num2str(tComp)]);

%verification
tic
violation30 = 0;
violation40 = 0;
violation50 = 0;
for iRun=1:length(Rcont)
    for i=1:length(Rcont{iRun})
        x_proj = interval(project(Rcont{iRun}{i},[1,4,7]));
        if any(infimum(x_proj) < -30)
            violation30 = 1;
        end
        if any(infimum(x_proj) < -40)
            violation40 = 1;
        end
        if any(infimum(x_proj) < -50)
            violation50 = 1;
        end
    end
end
violation30 
violation40 
violation50 
tVer = toc;
disp(['computation time of verification: ',num2str(tVer)]);


%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
totalInputChanges = 6;
fractionInputChange = 1;
tFinal = options.tFinal;

% set simulation options
stepsizeOptions = odeset('MaxStep',0.2*(options.tStart-options.tFinal));
% generate overall options
opt = odeset(stepsizeOptions);

% simulate all runs
for i=1:runs
    
    % initilaize results
    res.t{i} = [];
    res.x{i} = [];
    randInputCounter = 0;
    
    init = 1;
    
    %loop over discrete transitions
    for iTrans = 1:4
        %loop over input changes
        for iChange = 1:totalInputChanges
            %new run starts
            if init == 1
                %set initial state
                if i<=runs*fractionVertices
                    options.x0=randPointExtreme(Rinit); %initial state for simulation
                else
                    options.x0=randPoint(Rinit); %initial state for simulation
                end
                options.tStart = 0;
                options.tFinal = tFinal/totalInputChanges;
                init = 0;
            else
                options.tStart = options.tFinal;
                options.tFinal = options.tFinal + tFinal/totalInputChanges;
                options.x0 = x(end,:);
            end

            % set input
            %input from set of uncertain inputs 
            if randInputCounter <= fractionInputChange*iChange
                if i<=runs*fractionInputVertices
                    uRand = randPointExtreme(options.U); %random input from vertices
                else
                    uRand =randPoint(options.U); %random input from set
                end
                %update counter for changing random inputs
                randInputCounter = randInputCounter + 1;
            end
            % combine inputs
            options.u = uRand + options.uTrans; %input for simulation

            %simulate hybrid automaton
            if mod(iTrans,2)~=0
                [linSys_c,t,x] = simulate(linSys_c,options,options.tStart,options.tFinal,options.x0,opt); 
            else
                [linSys_n,t,x] = simulate(linSys_n,options,options.tStart,options.tFinal,options.x0,opt); 
            end
            res.t{i}(end+1:end+length(t),1) = t;
            res.x{i}(end+1:end+length(t),:) = x;
            if any(find(res.x{i} > 1e4))
                disp('stop');
            end
        end
    end
end


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
        for j=1:length(Rcont{i}) 
            Zproj = project(Rcont{i}{j},projectedDimensions);
            Zproj = reduce(Zproj,'girard',3);
            plotFilled(Zproj,[1 2],[.75 .75 .75],'EdgeColor','none');
        end
    end

    %plot initial set
    plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');
       
    %plot simulation results      
    for i=1:length(simRes.t)
        plot(res.x{i}(:,projectedDimensions(1)),res.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
    end

    %öabel plot
    xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
    ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
end

%plot results over time
for iDim = 1:1

    figure;
    hold on
    
    t1 = 0;
    t2 = options.timeStep;
    init = 1;

    %plot time elapse
    for i=1:length(Rcont)
        for j=1:length(Rcont{i}) 
            %get Uout 
            Uout1 = interval(project(Rcont{i}{j},iDim));
            %obtain times
            if init
                init = 0;
            else
                t1 = t1 + options.timeStep;
                t2 = t2 + options.timeStep;
            end
            %generate plot areas as interval hulls
            IH1 = interval([t1; infimum(Uout1)], [t2; supremum(Uout1)]);

            plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
        end
    end

    %plot simulation results
    for i=1:(length(res.t))
        plot(res.t{i},res.x{i}(:,iDim),'Color',0*[1 1 1]);
    end

end

%--------------------------------------------------------------------------

%------------- END OF CODE --------------