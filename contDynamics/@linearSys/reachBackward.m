function R = reachBackward(linsys,params,options)
% reachBackward - computes the AE or EA backward reachable set of a linear
%    continuous-time system, see [1]
%
% Syntax:
%    R = reachBackward(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    R - reachSet object storing the reachable set
%
% Example:
%    -
%
% References:
%    [1] M. Wetzlinger and M. Althoff, "Backward Reachability Analysis with
%        State Constraints for Linear Systems using Set Propagation", 2023.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:        Mark Wetzlinger
% Written:        04-December-2022
% Last update:    17-January-2023
% Last revision:  12-July-2023 (MW, move algs to private)

% ------------------------------ BEGIN CODE -------------------------------

% preprocessing
[params,options] = validateOptions(linsys,params,options);

% put sets in canonical form: U <- B*U
params.U = linsys.B * params.U + linsys.B * params.uTrans + linsys.c;

% disturbance already has to be of correct dimension (n)
params.W = linsys.E * params.W;

% first time-point is always the start set
timePoint.set{1} = params.R0;
timePoint.time{1} = 0;
timeInt = [];

% select correct algorithm
switch options.linAlg
    case {'inner:AE:timepoint','outer:AE:timepoint'}
        timePoint_ = reach_backward_AE_timepoint(linsys,params,options);
        timePoint.set = [timePoint.set; timePoint_.set];
        timePoint.time = [timePoint.time; timePoint_.time];
    case 'outer:AE:timeinterval'
        timeInt = reach_backward_AE_timeinterval(linsys,params,options);
    case {'inner:EA:timepoint','outer:EA:timepoint'}
        timePoint_ = reach_backward_EA_timepoint(linsys,params,options);
        timePoint.set = [timePoint.set; timePoint_.set];
        timePoint.time = [timePoint.time; timePoint_.time];
    case 'inner:EA:timeinterval'
        timeInt = reach_backward_EA_timeinterval(linsys,params,options);
end

% instantiate reachSet object
R = reachSet.initReachSet(timePoint,timeInt);

% ------------------------------ END OF CODE ------------------------------
