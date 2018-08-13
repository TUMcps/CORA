function res = test_nonlinear_reach_04_sevenDim_nonConvexRepr()
% test_nonlinear_reach_04_sevenDim_nonConvexRepr - unit_test_function of 
% nonlinear reachability analysis
%
% Checks the solution of a 7-dimensional nonlinear example using a non-convex
% set representation;
% It is checked whether the reachable set is enclosed in the initial set
% after a certain amount of time.
%
% Syntax:  
%    res = test_nonlinear_reach_04_sevenDim_nonConvexRepr()
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
% Written:      26-January-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


dim=7;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=0.2; %final time
%options.tFinal=6; %final time
options.x0=[1.2; 1.05; 1.5; 2.4; 1; 0.1; 0.45]; %initial state for simulation
options.R0=quadZonotope(options.x0,0.3*eye(dim),[],[],[]); %initial state for reachability analysis

options.timeStep=0.01; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=50; %zonotope order
options.intermediateOrder=5;
options.reductionTechnique='girard';
options.errorOrder=2;
options.polytopeOrder=10; %polytope order
options.reductionInterval=1e3;
options.maxError = 2*ones(dim,1);

options.plotType='frame';
options.projectedDimensions=[1 2];

options.originContained = 0;
options.advancedLinErrorComp = 1;
options.tensorOrder = 3;
%--------------------------------------------------------------------------


%obtain uncertain inputs
options.uTrans = 0;
options.U = zonotope([0]); %input for reachability analysis

%specify continuous dynamics-----------------------------------------------
sys=nonlinearSys(7,1,@sevenDimNonlinEq,options); %initialize system
%--------------------------------------------------------------------------


%compute reachable set using polynomial zonotopes
Rcont = reach(sys, options);

%enclose result by interval
IH = interval(Rcont{end}{1});

%saved result
IH_saved = interval( ...
    [1.0149039285033532;0.8007683032172868;0.9259718862023927;1.6088856146034343;0.4924670128231294;-0.0696710538369020;0.0041586523104461], ...
    [1.6994241453973689;1.5101149705032728;1.6826861040376888;2.4190315824628486;1.1030476700493148;0.2924851488683885;0.7099735189659464]);
      
        
%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_saved,1+1e-8));
res_2 = (IH_saved <= enlarge(IH,1+1e-8));

%final result
res = res_1*res_2;



%------------- END OF CODE --------------