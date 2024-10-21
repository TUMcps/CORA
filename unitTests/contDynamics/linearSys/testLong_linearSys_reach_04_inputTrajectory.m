function res = testLong_linearSys_reach_04_inputTrajectory
% testLong_linearSys_reach_04_inputTrajectory - unit test for linear
%    reachability  analysis with an input trajectory uTransVec; this test
%    should check whether correct the input trajectory is correctly
%    considered
%
% Syntax:
%    res = testLong_linearSys_reach_04_inputTrajectory
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       27-July-2018
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% load data
load data_Jean-MarcBiannic.mat A B uvec
dim_x = length(A);
dim_u = length(B(1,:));

% Parameters --------------------------------------------------------------

params.tFinal = 10; 
params.R0 = zonotope(zeros(dim_x,1),0.1*eye(dim_x,4));
params.u = [zeros(dim_u,1) uvec]; % input trajectory
params.u = params.u(:,1:1000);
% input for reachability analysis
params.U = zonotope(zeros(dim_u,1),diag([0.05,1]));


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.01;
options.taylorTerms = 4;
options.zonotopeOrder = 50;


% System Dynamics ---------------------------------------------------------

sys = linearSys('JeanMarcSys',A,B);


% Reachability Analysis (zonotope) ----------------------------------------

R = reach(sys, params, options);

% compute interval enclosure
IH = interval(R.timeInterval.set{end});

% saved result
IH_true = interval( ...
[-5.9848972828715841; -3.3983755171952517; -0.0000000000000028; -0.1350823568964674; -0.2148710415418167; -1.0946424167275999; -112.0033248270858337; -13.5774257689014330; 0.5841279576973016; 1.0309215393525788; -3.0015665115993233; -0.2755519893626190; 0.1114330836614218; -8.5036687640922963; -107.5757421162013259; -2.8820993171452698], ...
[5.5315325907849022; 3.8611903006005761; 0.0000000000000031; 0.2235395976356186; 0.6623604018319100; 0.7047701005021623; 117.3235376933028817; 1.3932307519464739; 0.7150847122290712; 1.2618079325715299; 1.7155254804071316; 0.4820624675390645; 0.6078848366564800; 11.9354856672451817; -7.0621927990805204; 0.7345517813137308]);
        
% check for equality up to tolerance
assert(isequal(IH,IH_true,1e-8));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
