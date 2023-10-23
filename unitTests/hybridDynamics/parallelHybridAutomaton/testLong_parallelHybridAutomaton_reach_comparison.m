function res = testLong_parallelHybridAutomaton_reach_comparison
% testLong_parallelHybridAutomaton_reach_comparison - unit test
%    function which compares the reachable set for a parallel hybrid
%    automaton with the one of the corresponding flat hybrid automaton
%    using the room heating benchmark described in Sec. 2.3 in [1] with two
%    rooms
%
% Syntax:
%    res = testLong_parallelHybridAutomaton_reach_comparison
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: roomHeating, roomHeatingParallel
%
% References: 
%   [1] A. Fehnker and F. Ivancic. "Benchmarks for Hybrid Systems 
%       Verification", HSCC 2004

% Authors:       Niklas Kochdumper
% Written:       29-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 3;
params.R0 = zonotope([20.5;20.5],diag([0.1,0.1]));  
params.U = zonotope(4,0.01);


% Reachability Settings ---------------------------------------------------

% settings for continuous reachability 
options.taylorTerms = 5; 
options.zonotopeOrder = 10; 
options.timeStep = 0.005; 
 
% settings for hybrid systems
options.enclose = {'box','pca'}; 
options.guardIntersect = 'zonoGirard';


% Hybrid Automaton --------------------------------------------------------

% system dynamics
HA = roomHeating();

% reachability analysis
params.startLoc = 1;
R_HA = reach(HA,params,options);


% Parallel Hybrid Automaton -----------------------------------------------

% system dynamics
PHA = roomHeatingParallel();

% reachability analysis
params.startLoc = [1;1];
R_PHA = reach(PHA,params,options);


% Numerical Evaluation ----------------------------------------------------

% interval enclosure of the final set
I_HA = interval(R_HA(end).timeInterval.set{end});
I_PHA = interval(R_PHA(end).timeInterval.set{end});

% check if slightly bloated versions enclose each other
res_1 = (I_HA <= enlarge(I_PHA,1+1e-10));
res_2 = (I_PHA <= enlarge(I_HA,1+1e-10));

% final result
res = res_1 && res_2;
    
% ------------------------------ END OF CODE ------------------------------
