function [params, R, simRes, unifiedOutputs, updatedSys] = confSynth(sys,params,options)
% confSynth - performs reachset conformance by computing the reachable 
%         of a system and checking whether all measurements are included
%
% Syntax:
%    [oarams, R, simRes, unifiedOutputs, updatedSys] = confSynth(sys,params,options)
%
% Inputs:
%    sys - contDynamics object
%    params - parameter defining the conformance problem
%    options - options for the conformance checking
%
% Outputs:
%    params - struct with conformant parameters
%    R - reachSet object (only time steps for which measurments exist)
%    simRes - states of the rapidly exploring random tree
%    unifiedOutputs - unified test cases (only deviation to nominal solution)
%    updatedSys - updated system
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% options preprocessing
options = validateOptions(sys, mfilename, params, options);

% initialize variables to be returned
R = [];
simRes = [];
unifiedOutputs = [];
updatedSys = [];

%% In case a specific algorithm working for any system dynamics has been chosen
% select conformance algorithm
if strcmp(options.confAlg,'RRT') % synthesize conformance using RRTs, see [1]
    [params, R, simRes] = confSynth_RRT(sys, params, options);
elseif strcmp(options.confAlg,'dyn') % conformance synthesis for dedicated system dynamics
    [params, unifiedOutputs] = confSynth_dyn(sys, params, options);
elseif strcmp(options.confAlg,'gray') % conformance synthesis for dedicated system dynamics
    [params, unifiedOutputs, updatedSys] = confSynth_gray(sys, params, options);
end

% ------------------------------ END OF CODE ------------------------------
