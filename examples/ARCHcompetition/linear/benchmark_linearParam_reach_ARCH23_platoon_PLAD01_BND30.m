function text = benchmark_linearParam_reach_ARCH23_platoon_PLAD01_BND30()
% benchmark_linearParam_reach_ARCH23_platoon_PLAD01_BND30 - example of linear
%    reachability analysis from the ARCH23 friendly competition
%    (platoon example); the linear dynamics switch deterministically
%
% Syntax:
%    benchmark_linearParam_reach_ARCH23_platoon_PLAD01_BND30
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
options.zonotopeOrder = 200;
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

% post-processing for plotting
t = 0;

for i = 1:size(R,1)
    % init time-interval struct
    timeInt.set = cell(length(R(i).timeInterval.set),1);
    timeInt.time = cell(length(R(i).timeInterval.set),1);

    for j = 1:length(R(i).timeInterval.set) 

        % compute interval
        timeInt.set{j} = interval(project(R(i).timeInterval.set{j},1));
        timeInt.time{j} = interval(t,t+options.timeStep);
        
        % update time
        t = t + options.timeStep;
    end

    % append to reachSet object
    R_new = reachSet([],timeInt);
    if i == 1
        R_full = R_new;
    else
        R_full = add(R_full,R_new);
    end
end


% Verification ------------------------------------------------------------

tic
violation30 = false;

for i = 1:size(R,1)
    for j = 1:length(R(i))
        x_proj = interval(project(R(i).timeInterval.set{j},[1,4,7]));
        if any(infimum(x_proj) < -30)
            violation30 = true;
        end
    end
end

tVer = toc;

disp(['specification verified: ',num2str(~violation30)]);
disp(['computation time: ',num2str(tVer+tComp)]);


% Simulation --------------------------------------------------------------

% parameters for simulation
runs = 10;
simOpt.points = 1;
simOpt.fracInpVert = 1;

% initial points
points = cell(runs,1);

% loop over all changes between controlled and uncontrolled system
x_all = cell(runs,1);
t_all = cell(runs,1);
y_all = cell(runs,1);

for j=1:runs

    points{j} = randPoint(R0);

    % reset time
    t = 0;
    simRes_ = {}; counter = 1;

    for i = 1:2
       
        % simulate controlled system
        params.tStart = t;
        params.tFinal = t + 5;
        t = t + 5;
        
        params.R0 = zonotope(points{j});
        temp = simulateRandom(linSys_c, params, simOpt);
        points{j} = temp.x{1}(end,:)';
        simRes_{counter} = temp;
        counter = counter + 1;
        
        % simulate uncontrolled system
        params.tStart = t;
        params.tFinal = t + 5;
        t = t + 5;
        
        params.R0 = zonotope(points{j});
        temp = simulateRandom(linSys_n, params, simOpt);
        points{j} = temp.x{1}(end,:)';
        simRes_{counter} = temp;
        counter = counter + 1;
    end

    % unify all parts to one simRes object
    x_all{j,1} = [simRes_{1}.x{1}; simRes_{2}.x{1}; simRes_{3}.x{1}; simRes_{4}.x{1}];
    t_all{j,1} = [simRes_{1}.t{1}; simRes_{2}.t{1}; simRes_{3}.t{1}; simRes_{4}.t{1}];
    y_all{j,1} = [simRes_{1}.y{1}; simRes_{2}.y{1}; simRes_{3}.y{1}; simRes_{4}.y{1}];

end

% init full simResult object
simRes = simResult(x_all,t_all,{},y_all);


% Visualization -----------------------------------------------------------

figure; hold on; box on;

useCORAcolors("CORA:contDynamics");

% plot reachable set
plotOverTime(R_full,1,'Unify',true);

updateColorIndex;

% plot simulation results
plotOverTime(simRes,1);

% labels
xlabel('$t$','interpreter','latex');
ylabel('$x_1$','interpreter','latex');


text = ['Platoon,PLAD01-BND30,',num2str(~violation30),',',num2str(tVer + tComp)];

% ------------------------------ END OF CODE ------------------------------
