function [res, R, simRes, unifiedOutputs] = confCheck(sys,params,options)
% confCheck - performs reachset conformance by computing the reachable 
%         of a system and checking whether all measurements are included
%
% Syntax:
%    [res, R, simRes, unifiedOutputs] = confCheck(sys,params,options)
%
% Inputs:
%    sys - contDynamics system
%    params - parameter defining the conformance problem
%    options - options for the conformance checking
%
% Outputs:
%    res - result: true/false 
%    R - reachSet object (only time steps for which measurments exist)
%    simRes - states of the rapidly exploring random tree
%    unifiedOutputs - unified test cases (only deviation to nominal solution)
%
% References:
%    [1] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.
%
% Example: 
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       15-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% options preprocessing
options = validateOptions(sys, mfilename, params, options);

% initialize variables to be returned
res = [];
R = [];
simRes = [];
unifiedOutputs = [];

%% In case a specific algorithm working for any system dynamics has been chosen
% select conformance algorithm
if strcmp(options.confAlg,'RRT') % check conformance using RRTs, see [1]
    [res, R, simRes] = confCheck_RRT(sys, params, options);
elseif strcmp(options.confAlg,'BF') % brute force conformance check (mainly for unit tests)
    [res, R] = confCheck_BF(sys, params, options);  
elseif strcmp(options.confAlg,'dyn') % conformance for dedicated system dynamics
    [res, unifiedOutputs] = confCheck_dyn(sys, params, options); 
end

% ------------------------------ END OF CODE ------------------------------
