function [res,R,simRes,unifiedOutputs] = isconform(sys,params,options,varargin)
% isconform - performs reachset conformance by computing the reachable set
%    of a system and checking whether all measurements are included
%
% Syntax:
%    [res,R,simRes,unifiedOutputs] = isconform(sys,params,options)
%    [res,R,simRes,unifiedOutputs] = isconform(sys,params,options,type)
%
% Inputs:
%    sys - contDynamics system
%    params - parameter defining the conformance problem
%    options - options for the conformance checking
%    type - "RRT", "BF" or "dyn", can be omitted for "dyn"
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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       15-June-2023
% Last update:   21-March-2024 (LL, change name from confCheck to isconform)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

type = setDefaultValues({"dyn"},varargin);

% options preprocessing
[params,options] = validateOptions(sys, params, options);

% initialize variables to be returned
R = [];
simRes = [];
unifiedOutputs = [];

%% In case a specific algorithm working for any system dynamics has been chosen
% select conformance algorithm
switch type
    case "RRT" % check conformance using RRTs, see [1]
        [res, R, simRes] = priv_isconform_RRT(sys, params, options);
    case "BF" % brute force conformance check (mainly for unit tests)
        [res, R] = priv_isconform_BF(sys, params, options);
    case "dyn" % conformance for dedicated system dynamics
        [res, unifiedOutputs] = isconform_dyn(sys, params, options);
    otherwise
        throw(CORAerror('CORA:wrongValue','fourth',...
            'type needs to be ''RRT'',''BF''  or ''dyn''.'))
end

% ------------------------------ END OF CODE ------------------------------
