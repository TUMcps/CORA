function res = test_requiredToolboxes
% test_requiredToolboxes - checks if required toolboxes are installed
%
% Syntax:
%    res = test_requiredToolboxes
%
% Inputs:
%    no
%
% Outputs:
%    res - true/false 
%
% Example: 
%

% Authors:       Matthias Althoff
% Written:       15-September-2016
% Last update:   04-May-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% note: license returns 1 if license exists and 0 if not (not: true/false!) 

% check if symbolic math toolbox is available
assert(logical(license('test','Symbolic_Toolbox')),'Symbolic Math Toolbox missing!');

% check if optimization toolbox is available
assert(logical(license('test','Optimization_Toolbox')),'Optimization Toolbox missing!');

% check if statistics toolbox is available
assert(logical(license('test','Statistics_Toolbox')),'Statistics Toolbox missing!');

p = path;

% note: contains returns true/false

% check if YALMIP toolbox is available
assert(isYalmipInstalled(),'YALMIP toolbox missing!');

% Check if at least one SDP solver supported by YALMIP is installed
assert(isSolverInstalled('mosek','sdpt3','gurobi','sedumi'),'SDP solver supported by YALMIP missing!');

% no check for MIP required as YALMIP 

% check if CORA is on the path
assert(contains(p,[filesep 'contDynamics' filesep]),'CORA is not on the Matlab path.'); 
assert(contains(p,[filesep 'contSet' filesep]),'CORA is not on the Matlab path.'); 
assert(contains(p,[filesep 'hybridDynamics' filesep]),'CORA is not on the Matlab path.');

res = true;

% ------------------------------ END OF CODE ------------------------------
