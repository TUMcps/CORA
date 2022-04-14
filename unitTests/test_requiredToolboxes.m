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
%    res - boolean 
%
% Example: 
%

% Author:       Matthias Althoff
% Written:      15-September-2016
% Last update:  04-May-2018
% Last revision:---

%------------- BEGIN CODE --------------

% note: license returns 1 if license exists and 0 if not (not: true/false!) 
% check if symbolic math toolbox is available
res_partial(1) = logical(license('test','Symbolic_Toolbox'));
if ~res_partial(1)
    disp('Symbolic Toolbox missing!');
end

% check if optimization toolbox is available
res_partial(2) = logical(license('test','Optimization_Toolbox'));
if ~res_partial(2)
    disp('Optimization Toolbox missing!');
end

res_partial(3) = logical(license('test','Statistics_Toolbox'));
if ~res_partial(3)
    disp('Statistics Toolbox missing!');
end

p = path;

% note: contains returns true/false
% check if MPT toolbox is available
res_partial(4) = contains(p,[filesep 'mpt' filesep]);
if ~res_partial(4)
    disp('MPT toolbox missing!');
end

% check if YALMIP toolbox is available
res_partial(5) = contains(p,[filesep 'yalmip' filesep]);
if ~res_partial(5)
    disp('YALMIP toolbox missing!');
end

% check if CORA is on the path
res_partial(6) = contains(p,[filesep 'contDynamics' filesep]); 
res_partial(7) = contains(p,[filesep 'contSet' filesep]); 
res_partial(8) = contains(p,[filesep 'hybridDynamics' filesep]);

res = all(res_partial);

%------------- END OF CODE --------------
