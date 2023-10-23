function res = test_linParamSys_reach_02_5dim_var()
% test_linParamSys_reach_02_5dim_var - unit test of linear parametric
%    reachability analysis from [1] where the parameters vary over time
%
% Syntax:
%    res = test_linParamSys_reach_02_5dim_var
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] Althoff, M.; Le Guernic, C. & Krogh, B. H. Reachable Set Computation
%        for Uncertain Time-Varying Linear Systems. HSCC 2011, 93-102

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       03-October-2017
% Last update:   23-April-2020 (restructure params/options)
%                05-June-2020 (NK, adapted to bug fix in commit fdc7bba)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameters --------------------------------------------------------------

dim_x = 5;
params.R0 = zonotope(ones(dim_x,1),0.1*eye(dim_x)); %initial state for reachability analysis
params.U = zonotope(zeros(dim_x,1),0.1*eye(dim_x)); %input for reachability analysis
params.tFinal = 5; %final time


% Reachability Settings ---------------------------------------------------

% (taken from reference)
options.timeStep = 0.05; %time step size for reachable set computation
options.taylorTerms = 4; %number of taylor terms for reachable sets
options.intermediateTerms = 2;
options.zonotopeOrder = 20; %zonotope order


% System Dynamics ---------------------------------------------------------

Acenter = [-1 -4  0  0  0;
            4 -1  0  0  0;
            0  0 -3  1  0;
            0  0 -1 -3  0;
            0  0  0  0 -2];
Arad{1} = [0.1 0.1 0   0   0;
           0.1 0.1 0   0   0; 
           0   0   0.1 0.1 0;
           0   0   0.1 0.1 0;
           0   0   0   0   0.1];
matZ_A = matZonotope(Acenter,Arad);
matI_A = intervalMatrix(matZ_A);

fiveDimSys_zono = linParamSys(matZ_A, 1,'varParam'); %instantiate system
fiveDimSys_int = linParamSys(matI_A, 1,'varParam'); %instantiate system


% Reachability Analysis ---------------------------------------------------

% compute reachable set using matrix zonotope / interval matrix
R_zono = reach(fiveDimSys_zono, params, options);
R_int = reach(fiveDimSys_int, params, options);


% Numerical Evaluation ----------------------------------------------------

IH_zono = interval(R_zono.timeInterval.set{end});
IH_int = interval(R_int.timeInterval.set{end});

%saved result
IH_saved_zono = interval( ...
        [ -0.211200438431055; -0.200109643200870; -0.052530295054846; -0.052626609458356; -0.058823653515915], ...
        [0.206114142337843; 0.218943349986004; 0.052529828329088; 0.052627418115272; 0.058919155036720]);
IH_saved_int = interval( ...
        [  -0.258708023885910; -0.246551871520107; -0.054033103506035; -0.054079269389896; -0.058822997090575], ...
        [0.253670265974255; 0.265394513562644; 0.054032634687992; 0.054080076731048; 0.058918625877829]);

%final result
res = isequal(IH_zono,IH_saved_zono,1e-8) && isequal(IH_int,IH_saved_int,1e-8);

% ------------------------------ END OF CODE ------------------------------
