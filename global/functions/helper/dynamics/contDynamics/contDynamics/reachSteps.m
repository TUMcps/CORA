function val = reachSteps(params,options)
% reachSteps - computes number of steps in execution of reachability
%    algorithm for a given time horizon and time step size
%
% Syntax:
%    val = reachSteps(params,options)
%
% Inputs:
%    params - model parameters
%    options - reachability settings
%
% Outputs:
%    val - number of steps in execution of reachability algorithm
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Victor Gassmann
% Written:       15-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

val = round((params.tFinal-params.tStart)/options.timeStep);

% ------------------------------ END OF CODE ------------------------------
