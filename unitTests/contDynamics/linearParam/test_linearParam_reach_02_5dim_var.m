function res = test_linearParam_reach_02_5dim_var()
% test_linearParam_reach_02_5dim_var - unit test of 
% linear parametric reachability analysis where the parameters can vary 
% over time.
%
% This example is taken from
%
% Althoff, M.; Le Guernic, C. & Krogh, B. H. Reachable Set Computation for 
% Uncertain Time-Varying Linear Systems Hybrid Systems: Computation and 
% Control, 2011, 93-102
%
% Syntax:  
%    res = test_linearParam_reach_02_5dim_var
%
% Inputs:
%    none
%
% Outputs:
%    none
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      03-October-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


dim = 5;

%set options --------------------------------------------------------------
options.tStart = 0; %start time
options.tFinal = 5; %final time
options.x0 = ones(dim,1); %initial state for simulation
options.R0 = zonotope([options.x0,0.1*eye(length(options.x0))]); %initial state for reachability analysis

options.timeStep = 0.05; %time step size for reachable set computation
options.taylorTerms = 4; %number of taylor terms for reachable sets
options.intermediateOrder = 2;
options.zonotopeOrder = 20; %zonotope order
options.originContained = 0;
options.reductionTechnique = 'girard';

uTrans = zeros(dim,1);
options.uTrans = uTrans; %center of input set
options.U = zonotope([zeros(dim,1),0.1*eye(dim)]); %input for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
Acenter = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
Arad{1} = [0.1 0.1 0 0 0; 0.1 0.1 0 0 0; 0 0 0.1 0.1 0; 0 0 0.1 0.1 0; 0 0 0 0 0.1];
matZ_A = matZonotope(Acenter,Arad);
matI_A = intervalMatrix(matZ_A);

fiveDimSys_zono = linParamSys(matZ_A, 1, options.timeStep, options.taylorTerms,'varParam'); %instantiate system
fiveDimSys_int = linParamSys(matI_A, 1, options.timeStep, options.taylorTerms,'varParam'); %instantiate system
%--------------------------------------------------------------------------

%compute reachable set using zonotopes
Rcont_zono = reach(fiveDimSys_zono, options);
Rcont_int = reach(fiveDimSys_int, options);

IH_zono = interval(Rcont_zono{end});
IH_int = interval(Rcont_int{end});

%saved result
IH_saved_zono = interval( ...
        [-0.210884705282954; -0.200286600390145; -0.052519883462392; -0.052623375987781; -0.058823568694747], ...
        [0.205837073920494; 0.219081898999744; 0.052519415273674; 0.052624185268516; 0.058919135578997]);
IH_saved_int = interval( ...
        [-0.258418678296276; -0.246526512465825; -0.054008238679980; -0.054071241145380; -0.058822912184395], ...
        [0.253420403100503; 0.265329676949767; 0.054007768375767; 0.054072049089100; 0.058918606421781]);

%check if slightly bloated versions enclose each other for IH
res_11 = (IH_zono <= enlarge(IH_saved_zono, 1+1e-8));
res_12 = (IH_saved_zono <= enlarge(IH_zono, 1+1e-8));

%check if slightly bloated versions enclose each other for IH2
res_21 = (IH_int <= enlarge(IH_saved_int, 1+1e-8));
res_22 = (IH_saved_int <= enlarge(IH_int, 1+1e-8));

%final result
res = res_11*res_12*res_21*res_22;

%------------- END OF CODE --------------