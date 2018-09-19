function res = test_linear_reach_05_inputTrajectoryOnly()
% test_linear_reach_05_inputTrajectoryOnly - unit test for linear 
% reachability analysis with an input trajectory uTransVec; this test 
% should check whether correct the input trajectory is correctly considered
%
% Syntax:  
%    res = test_linear_reach_05_inputTrajectoryOnly
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
% Written:      18-September-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

dim=5;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=0.12; %final time
options.x0=ones(dim,1); %initial state for simulation
options.R0=zonotope([options.x0,0*eye(length(options.x0))]); %initial state for reachability analysis

options.timeStep=0.04; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=200; %zonotope order
options.originContained=0;
options.reductionTechnique='girard';

options.uTransVec = 100*[1; 0; 0; 0.5; -0.5]; % start of input trajectory
options.uTransVec(:,2:3) = 0;
options.U=0*zonotope([zeros(5,1),diag([0.2, 0.5, 0.2, 0.5, 0.5])]); %input for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
A=[-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B=1;
fiveDimSys=linearSys('fiveDimSys',A,B); %initialize system
%--------------------------------------------------------------------------

% compute reachable set using zonotopes
Rcont = reach(fiveDimSys, options);

% interval hull of final set
IH = interval(Rcont{end});

%saved result
IH_true = interval( ...
[3.7035909068940835; 2.0587487720830868; 0.9218120355943664; 2.0792489149905693; -0.9222005200165104], ...
[4.2545447422606566; 2.6163966074240355; 0.9494850143887323; 2.3899927521345345; -0.8505317877218803]);
        
%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_true,1+1e-8));
res_2 = (IH_true <= enlarge(IH,1+1e-8));

%final result
res = res_1*res_2;

%------------- END OF CODE --------------
