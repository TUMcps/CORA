function example_linearParam_reach_ARCH21_platoon_PLAA01_BND42()
% example_linearParam_reach_ARCH21_platoon_PLAA01_BND42 - example of linear 
%    reachability analysis from the ARCH21 friendly competition
%    (platoon example); the linear dynamics can switch arbitrarily
%
% Syntax:  
%    example_linearParam_reach_ARCH21_platoon_PLAA01_BND42
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 

% Author:       Matthias Althoff
% Written:      09-February-2017
% Last update:  13-March-2019
% Last revision:---

%------------- BEGIN CODE --------------


% Parameter ---------------------------------------------------------------

% initial set
params.R0 = zonotope(zeros(9,1));

% uncertain inputs
params.U = zonotope(interval(-9,1));

% final time
params.tFinal = 20; 


% Reachability Settings ---------------------------------------------------

options.intermediateTerms = 3;
options.timeStep = params.tFinal/2222;
options.zonotopeOrder = 800;
options.taylorTerms = 3;
options.compTimePoint = false;


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
        
% build zonotpe matrix
A_mid = 0.5*(A_c + A_n);
A_gen{1} = 0.5*(A_c - A_n);
matZ_A = matZonotope(A_mid, A_gen);

% instantiate linear dynamics with constant parameters
linSys  = linParamSys(matZ_A, B,'constParam');


% Reachability Analysis ---------------------------------------------------

tic
R = reach(linSys, params, options);
tComp = toc;



% Verification ------------------------------------------------------------

tic
violation42 = 0;
for i=1:length(R.timeInterval.set)
    x_proj = interval(project(R.timeInterval.set{i},[1,4,7]));
    if any(infimum(x_proj) < -42)
        violation42 = 1;
    end
end
tVer = toc;

disp(['specifications verified: ',num2str(~violation42)]);
disp(['computation time: ',num2str(tVer+tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 60;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;

simRes = simulateRandom(linSys, params, simOpt);



% Visualization -----------------------------------------------------------

figure; hold on; box on;

% plot reachable set over time 
plotOverTime(R,1,'FaceColor',[.75 .75 .75],'EdgeColor','none');

% plot simulation results
plotOverTime(simRes,1);

% formatting
xlabel('t');
ylabel('x_1');

%------------- END OF CODE --------------