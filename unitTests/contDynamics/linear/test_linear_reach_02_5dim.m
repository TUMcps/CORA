function res = test_linear_reach_02_5dim()
% test_linear_reach_02_5dim - unit_test_function of linear reachability 
% analysis with uncertain inputs
%
% Checks the solution of the linearSys class for a 5-dimensional example;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:  
%    res = test_linear_reach_02_5dim()
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
% Written:      09-August-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

dim=5;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=5; %final time
options.x0=ones(dim,1); %initial state for simulation
options.R0=zonotope([options.x0,0.1*eye(length(options.x0))]); %initial state for reachability analysis

options.timeStep=0.04; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=10; %zonotope order
options.originContained=0;
options.reductionTechnique='girard';
options.linAlg = 1;

uTrans=[1; 0; 0; 0.5; -0.5];
options.uTrans=uTrans; %input for simulation
options.U=0.5*zonotope([zeros(5,1),diag([0.2, 0.5, 0.2, 0.5, 0.5])]); %input for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
A=[-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B=1;
fiveDimSys=linearSys('fiveDimSys',A,B); %initialize system
%--------------------------------------------------------------------------


%compute reachable set using zonotopes
Rcont = reach(fiveDimSys, options);

IH = interval(Rcont{end});

%saved result
IH_saved = interval( ...
           [-0.186078149456309; -0.004994826957236; -0.010811262167507; 0.053366186432215; -0.385353029993981], ...
           [0.300725848257715; 0.491957694383420; 0.110810877781333; 0.246634559097104; -0.114528793200091]);
        
%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_saved,1+1e-8));
res_2 = (IH_saved <= enlarge(IH,1+1e-8));

%final result
res_zono = res_1*res_2;


%use zonotope bundles instead of zonotopes---------------------------------
%generate zonotope bundle
Z0{1} = options.R0;
Z0{2} = options.R0 + [-0.1; 0; 0.1; 0; 0];
options.R0 = zonotopeBundle(Z0);

%compute reachable set using zonotope bundles
Rcont = reach(fiveDimSys, options);

IH = interval(Rcont{end});

%saved result
IH_saved = interval( ...
    [-0.1860781494563091; -0.0049948269572356; -0.0108112537143071; 0.0533662157659313; -0.3853530299939805], ...
    [0.3003413164189264; 0.4913712245778811; 0.1108108777813331; 0.2466345590971040; -0.1145287932000914]);
        
%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_saved,1+1e-8));
res_2 = (IH_saved <= enlarge(IH,1+1e-8));

%final result
res_zonoBundles = res_1*res_2;

%result of different set representations
res = res_zono*res_zonoBundles;

%------------- END OF CODE --------------
