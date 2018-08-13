function example_linear_reach_ARCH18_building()
% example_linear_reach_ARCH18_building - example of linear reachability 
% analysis from the ARCH18 friendly competition (building example with time
% varying inputs)
%
% Syntax:  
%    example_linear_reach_ARCH18_building
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
% Written:      28-May-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

load build
dim=length(A);



% Options -----------------------------------------------------------------

R0 = interval([0.0002*ones(10,1); -0.0001; zeros(37,1)], [0.00025*ones(10,1); 0.0001; zeros(37,1)]);
options.x0=mid(R0); %initial state for simulation

options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=100; %zonotope order
options.originContained=0;
options.reductionTechnique='girard';
options.linAlg = 1;

uTrans=0.9;
options.uTrans=uTrans; %center of input set
options.U=zonotope([options.uTrans,0.1]); %input for reachability analysis

buildingSys=linearSys('buildingSys',A,B); %initialize system




% Reachability Analysis ---------------------------------------------------

%compute reachable set using zonotopes
tic
%1st phase
options.R0=zonotope(R0); %initial state for reachability analysis
options.tStart=0; %start time
options.tFinal=1; %final time
timeStep_1 = 0.002;
options.timeStep = timeStep_1;
Rcont = reach(buildingSys, options);
%2ndphase
options.R0 = Rcont{end};
options.tStart=1; %start time
options.tFinal=20; %final time
timeStep_2 = 0.01;
options.timeStep = timeStep_2;
Rcont2 = reach(buildingSys, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);




% Verification ------------------------------------------------------------

tic
violation = 0;
boundReached1 = 0;
boundReached2 = 0;

for i=1:length(Rcont)
    x_25 = interval(project(Rcont{i},25));
    if supremum(x_25) > 5.1e-3
        violation = 1;
    end
    if supremum(x_25) > 4e-3
        boundReached1 = 1;
    end
end
for i=1:length(Rcont2)
    x_25 = interval(project(Rcont2{i},25));
    if supremum(x_25) > 5.1e-3
        violation = 1;
    end
end
x_25 = interval(project(Rcont2{end},25));
if infimum(x_25) < -0.78e-3
   boundReached2 = 1; 
end

violation
boundReached1
boundReached2

tVer = toc;
disp(['computation time of verification: ',num2str(tVer)]);




% Simulation --------------------------------------------------------------

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 10;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 20;
options.R0=zonotope(R0); %initial state for simulation results
simRes = simulate_random(buildingSys, options, runs, fractionVertices, fractionInputVertices, inputChanges);





% Visualization -----------------------------------------------------------

% Plot 1: time interval t \in [0,1] s

figure;
hold on

i = 1;
iDim = 25;
%plot time elapse
while i<=length(Rcont)
    %get Uout 
    t1 = (i-1)*timeStep_1;
    try
        minVal = inf;
        maxVal = -inf;
        for k=1:1
            Uout = interval(project(Rcont{i},iDim));
            if infimum(Uout) < minVal
                minVal = infimum(Uout);
            end
            if supremum(Uout) > maxVal
                maxVal = supremum(Uout);
            end
            i = i + 1;
        end
    catch
        minVal = infimum(interval(project(Rcont{i-1},iDim)));
        maxVal = supremum(interval(project(Rcont{i-1},iDim)));
    end
    t2 = (i-1)*timeStep_1;
    %generate plot areas as interval hulls
    IH1 = interval([t1; minVal], [t2; maxVal]);

    plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
end

% Simulation
for i=1:(length(simRes.t))
    plot(simRes.t{i},simRes.x{i}(:,25),'Color',0*[1 1 1]);
end

axis([0, 1, -6.5e-3, 6e-3])
xlabel('t');
ylabel('x_{25}');
box on



% Plot 2: time interval t \in [0,20] s

figure;
hold on

i = 1;
iDim = 25;

% First part 
while i<=length(Rcont)
    %get Uout 
    t1 = (i-1)*timeStep_1;
    try
        minVal = inf;
        maxVal = -inf;
        for k=1:10
            Uout = interval(project(Rcont{i},iDim));
            if infimum(Uout) < minVal
                minVal = infimum(Uout);
            end
            if supremum(Uout) > maxVal
                maxVal = supremum(Uout);
            end
            i = i + 1;
        end
    catch
        minVal = infimum(interval(project(Rcont{i-1},iDim)));
        maxVal = supremum(interval(project(Rcont{i-1},iDim)));
    end
    t2 = (i-1)*timeStep_1;
    %generate plot areas as interval hulls
    IH1 = interval([t1; minVal], [t2; maxVal]);

    plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
end

% Second part
i = 1;
while i<=length(Rcont2)
    %get Uout 
    t1 = options.tStart + (i-1)*timeStep_2;
    try
        minVal = inf;
        maxVal = -inf;
        for k=1:5
            Uout = interval(project(Rcont2{i},iDim));
            if infimum(Uout) < minVal
                minVal = infimum(Uout);
            end
            if supremum(Uout) > maxVal
                maxVal = supremum(Uout);
            end
            i = i + 1;
        end
    catch
        minVal = infimum(interval(project(Rcont2{i-1},iDim)));
        maxVal = supremum(interval(project(Rcont2{i-1},iDim)));
    end
    t2 = options.tStart + (i-1)*timeStep_2;
    %generate plot areas as interval hulls
    IH1 = interval([t1; minVal], [t2; maxVal]);

    plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
end

% Simulation
for i=1:(length(simRes.t))
    plot(simRes.t{i},simRes.x{i}(:,25),'Color',0*[1 1 1]);
end

axis([0, 20, -6.5e-3, 6e-3])
xlabel('t');
ylabel('x_{25}');
box on


%------------- END OF CODE --------------
