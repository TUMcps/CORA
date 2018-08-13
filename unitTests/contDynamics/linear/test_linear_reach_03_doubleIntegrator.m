function res = test_linear_reach_03_doubleIntegrator()
% test_linear_reach_03_doubleIntegrator - unit test for linear reachability 
% analysis with uncertain inputs; this test should check whether correct
% results are returned when the system matrix only consists of zeros
%
% Syntax:  
%    res = test_linear_reach_03_doubleIntegrator()
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
% Written:      12-November-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

dim = 2;

%set options --------------------------------------------------------------
options.tStart = 0; %start time
options.tFinal = 1; %final time
options.x0 = ones(dim,1); %initial state for simulation
options.R0 = zonotope([options.x0,0.1*eye(length(options.x0))]); %initial state for reachability analysis

options.timeStep = 0.04; %time step size for reachable set computation
options.taylorTerms = 4; %number of taylor terms for reachable sets
options.zonotopeOrder = 10; %zonotope order
options.originContained = 0;
options.reductionTechnique = 'girard';
options.linAlg = 1;

uTrans = [1; 1];
options.uTrans = uTrans; %input for simulation
options.U = zonotope([zeros(2,1),diag([0.1, 0.1])]); %input for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
A = zeros(2);
B = 1;
doubleIntegrator = linearSys('twoDimSys',A,B); %initialize system
%--------------------------------------------------------------------------

%compute reachable set using zonotopes
Rcont = reach(doubleIntegrator, options);

IH = interval(Rcont{end});

%saved result
IH_true = ones(dim,1)*interval(1.76, 2.2);
        
%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_true,1+1e-8));
res_2 = (IH_true <= enlarge(IH,1+1e-8));

%final result
res = res_1*res_2;

%------------- END OF CODE --------------
