function text = benchmark_linearParam_reach_ARCH23_platoon_PLAD01_BND42()
% benchmark_linearParam_reach_ARCH23_platoon_PLAD01_BND42 - example of linear
%    reachability analysis from the ARCH23 friendly competition
%    (platoon example); the linear dynamics switch deterministically
%
% Syntax:
%    benchmark_linearParam_reach_ARCH23_platoon_PLAD01_BND42
%
% Inputs:
%    -
%
% Outputs:
%    -

% Authors:       Matthias Althoff
% Written:       05-April-2017
% Last update:   13-March-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameter ---------------------------------------------------------------

% initial set
R0 = zonotope(zeros(9,1));

% uncertain inputs
params.U = zonotope(interval(-9,1));

% final time
params.tFinal = 5;


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.02; 
options.zonotopeOrder = 20;
options.taylorTerms = 4;


% System Dynamics ---------------------------------------------------------

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

% instantiate linear dynamics with constant parameters
linSys_c  = linearSys('c',A_c, B);
linSys_n  = linearSys('n',A_n, B);


% Reachability Analysis ---------------------------------------------------

clock = tic;
R = []; 
params.R0 = R0;

% loop over all changes between controlled and uncontrolled system
for i = 1:2    
    % reachable set for controlled system
    Rtemp = reach(linSys_c, params, options);
    
    params.R0 = Rtemp.timePoint.set{end};
    R = add(R,Rtemp);
    
    % reachable set for uncnotrolled system
    Rtemp = reach(linSys_n, params, options);
    params.R0 = Rtemp.timePoint.set{end};
    R = add(R,Rtemp);
end

tComp = toc(clock);


% Verification ------------------------------------------------------------

clock = tic;
violation42 = false;

for i = 1:size(R,1)
    for j = 1:length(R(i).timeInterval.set)
        x_proj = interval(project(R(i).timeInterval.set{j},[1,4,7]));
        if any(infimum(x_proj) < -42)
            violation42 = true;
        end
    end
end

tVer = toc(clock);

disp(['specification verified: ',num2str(~violation42)]);
disp(['computation time: ',num2str(tVer+tComp)]);


% Simulation --------------------------------------------------------------

% parameter for simulation
runs = 60;

simOpt.points = 1;
simOpt.fracInpVert = 1;

% initial points
points = cell(runs,1);
simRes = cell(4*runs,1);

for i = 1:runs
    points{i} = randPoint(R0); 
end

% loop over all changes between controlled and uncontrolled system
t = 0;
counter = 1;

for i = 1:2
   
    % simulate controlled system
    if t ~= 0
        params.tStart = t;
    end
    params.tFinal = t + 5;
    t = t + 5;
    
    for j = 1:length(points)
        params.R0 = zonotope(points{j});
        temp = simulateRandom(linSys_c, params, simOpt);
        points{j} = temp.x{1}(end,:)';
        simRes{counter} = temp;
        counter = counter + 1;
    end
    
    % simulate uncontrolled system
    params.tStart = t;
    params.tFinal = t + 5;
    t = t + 5;
    
    for j = 1:length(points)
        params.R0 = zonotope(points{j});
        temp = simulateRandom(linSys_n, params, simOpt);
        points{j} = temp.x{1}(end,:)';
        simRes{counter} = temp;
        counter = counter + 1;
    end
end


% Visualization -----------------------------------------------------------

figure; hold on, box on;

% reachable set over time
t = 0;

for i = 1:size(R,1)
    for j = 1:length(R(i).timeInterval.set) 

        % compute interval
        intX = interval(project(R(i).timeInterval.set{j},1));
        intT = interval(t,t+options.timeStep);
        int = cartProd(intT,intX);
        
        % plot interval
        plot(int,[1 2],'FaceColor',CORAcolor("CORA:reachSet"));
        
        % update time
        t = t + options.timeStep;
    end
end

% plot simulation results
for i = 1:length(simRes)
    for j = 1:length(simRes{i}.t)
        plot(simRes{i}.t{j},simRes{i}.x{j}(:,1),'Color',CORAcolor("CORA:simulations"));
    end
end

xlabel('t');
ylabel('x_1');


text = ['Platoon,PLAD01-BND42,',num2str(~violation42),',',num2str(tVer + tComp)];

% ------------------------------ END OF CODE ------------------------------
