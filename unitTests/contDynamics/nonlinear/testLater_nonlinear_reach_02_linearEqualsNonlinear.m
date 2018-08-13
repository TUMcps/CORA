function res = test_nonlinear_reach_02_linearEqualsNonlinear()
% test_nonlinear_reach_02_linearEqualsNonlinear - unit_test_function for 
% nonlinear reachability analysis; it is checked whether the solution of a
% linear system equals the solution obtained by the nonlinear sytsem class
% when the system is linear
%
% Compares the solution of the linearSys class for a 5-dimensional linear 
% example; It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:  
%    res = test_nonlinear_reach_02_linearEqualsNonlinear()
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
options.taylorTerms=8; %number of taylor terms for reachable sets
options.zonotopeOrder=10; %zonotope order
options.intermediateOrder=options.zonotopeOrder;
options.originContained=0;
options.reductionTechnique='girard';
options.errorOrder=1;
options.tensorOrder = 2;
options.reductionInterval=1e3;
options.maxError = 1*ones(dim,1);
options.linAlg = 1;

uTrans=[1; 0; 0; 0.5; -0.5];
options.uTrans=uTrans; %input for simulation
options.U=0.5*zonotope([zeros(5,1),diag([0.2, 0.5, 0.2, 0.5, 0.5])]); %input for reachability analysis
%--------------------------------------------------------------------------

%obtain continuous linear dynamics-----------------------------------------
A=[-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B=1;
fiveDimSys=linearSys('fiveDimSys',A,B); %initialize system
%--------------------------------------------------------------------------

%write as equivalent nonlinear dynamics------------------------------------
fiveDimSysNonlinear=nonlinearSys(5,5,@fiveDimSysEq,options); %initialize system
%--------------------------------------------------------------------------

%compute reachable set using linear dynamics class
Rcont = reach(fiveDimSys, options);
%compute reachable set using nonlinear dynamics class (tensor order 1)
options.advancedLinErrorComp = 0;
options.tensorOrder = 1;
RcontNonlinear_T1 = reach(fiveDimSysNonlinear, options);
%compute reachable set using nonlinear dynamics class (tensor order 2)
options.advancedLinErrorComp = 1;
options.tensorOrder = 2;
RcontNonlinear_T2 = reach(fiveDimSysNonlinear, options);

%enclosing intervals of final reachable sets
IH = interval(Rcont{end});
IH_nonlinear_T1 = interval(RcontNonlinear_T1{end}{1});
IH_nonlinear_T2 = interval(RcontNonlinear_T2{end}{1});
        
%check if slightly bloated versions enclose each other (linear vs T1)
res_T1_1 = (IH <= enlarge(IH_nonlinear_T1,1+1e-8));
res_T1_2 = (IH_nonlinear_T1 <= enlarge(IH,1+1e-8));

%check if slightly bloated versions enclose each other (linear vs T2)
res_T2_1 = (IH <= enlarge(IH_nonlinear_T2,1+1e-8));
res_T2_2 = (IH_nonlinear_T2 <= enlarge(IH,1+1e-8));

%final result
res = res_T1_1*res_T1_2*res_T2_1*res_T2_2;

%------------- END OF CODE --------------
