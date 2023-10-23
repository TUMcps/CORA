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

Res = [];
% note: license returns 1 if license exists and 0 if not (not: true/false!) 
% check if symbolic math toolbox is available
Res = logical(license('test','Symbolic_Toolbox'));
if ~Res(end)
    disp('Symbolic Math Toolbox missing!');
end

% check if optimization toolbox is available
Res = [Res,logical(license('test','Optimization_Toolbox'))];
if ~Res(end)
    disp('Optimization Toolbox missing!');
end

Res = [Res,logical(license('test','Statistics_Toolbox'))];
if ~Res(end)
    disp('Statistics Toolbox missing!');
end

p = path;

% note: contains returns true/false

% check if YALMIP toolbox is available
Res = [Res,isYalmipInstalled()];
if ~Res(end)
    disp('YALMIP toolbox missing!');
end

% Check if at least one SDP solver supported by YALMIP is installed
Res = [Res,isSolverInstalled('mosek','sdpt3','gurobi','sedumi')];
if ~Res(end)
    disp('SDP solver supported by YALMIP missing!');
end

% no check for MIP required as YALMIP 


% check if CORA is on the path
Res = [Res,contains(p,[filesep 'contDynamics' filesep])]; 
Res = [Res,contains(p,[filesep 'contSet' filesep])]; 
Res = [Res,contains(p,[filesep 'hybridDynamics' filesep])];

res = all(Res);

% ------------------------------ END OF CODE ------------------------------
