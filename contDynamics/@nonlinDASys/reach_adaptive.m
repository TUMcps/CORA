function [timeInt,timePoint,res,tVec,options] = reach_adaptive(nlnsysDA,params,options)
% reach_adaptive - computes the reachable set for a nonlinear
%    differential-algebraic system using adaptive parameter tuning
%
% Syntax:
%    [timeInt,timePoint,res,tVec,options] = reach_adaptive(nlnsysDA,params,options)
%
% Inputs:
%    nlnsysDA - nonlinDASys object
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    timeInt - cell-array of time-interval solutions
%    timePoint - cell-array of time-point solutions
%    res - satisfaction / violation of specifications
%    tVec - vector of time steps
%    options - options for the computation of reachable sets (param tracking)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       17-June-2021
% Last update:   30-August-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize cell-arrays that store the reachable set
timeInt.set = {};
timeInt.time = {};
timeInt.algebraic = {};
timePoint.set = {};
timePoint.time = {};
res = true;
tVec = 0;

% remove 'adaptive' from alg (just for tensor computation)
if contains(options.alg,'lin')
    options.alg = 'lin';
end

% iteration counter and time for main loop
options.i = 1;
options.t = params.tStart;
options.R = params.R0;

% obtain consistent initial algebraic set
y0 = params.y0guess;
y0 = consistentInitialState(nlnsysDA, center(options.R), y0, params.uTrans);
Rstart_y = zonotope(y0);
% set linearization errors (separate values for Delta and optimal Delta t)
options.error_adm_x_horizon = eps*ones(nlnsysDA.nrOfStates,1);
options.error_adm_x_Deltatopt = eps*ones(nlnsysDA.nrOfStates,1);
options.error_adm_y_horizon = eps*ones(nlnsysDA.nrOfConstraints,1);
options.error_adm_y_Deltatopt = eps*ones(nlnsysDA.nrOfConstraints,1);
% init abortion flag
abortAnalysis = false;

% MAIN LOOP
while params.tFinal - options.t > 1e-12 && ~abortAnalysis
    
    % log information
    verboseLog(options.verbose,options.i,options.t,params.tStart,params.tFinal);
    
    % reachable set propagation
    [Rnext.ti,Rnext.tp,Rstart_y,options] = linReach_adaptive(nlnsysDA,options.R,Rstart_y,params,options);
    
    % reduction for next step
    Rnext.ti = reduce(Rnext.ti,'adaptive',options.redFactor*5); % not reused
    Rnext.tp = reduce(Rnext.tp,'adaptive',options.redFactor);
    

    % save to output variables
    tVec(options.i,1) = options.timeStep;
    % save reachable set in cell structure
    timeInt.set{options.i,1} = Rnext.ti; 
    timeInt.time{options.i,1} = interval(options.t,options.t+tVec(options.i));
    timeInt.algebraic{options.i,1}  = Rstart_y;
    timePoint.set{options.i,1} = Rnext.tp;
    timePoint.time{options.i,1} = options.t+tVec(options.i);
    
    % increment time
    options.t = options.t + options.timeStep;
    
    % update iteration counter
    options.i = options.i + 1;
    
    % start set for next step (since always initReach called)
    options.R = Rnext.tp;
    
    % check for timeStep -> 0
    abortAnalysis = aux_checkForAbortion(tVec,options.t,params.tFinal);
end

% log information
verboseLog(options.verbose,options.i,options.t,params.tStart,params.tFinal);

end


% Auxiliary functions -----------------------------------------------------

function abortAnalysis = aux_checkForAbortion(tVec,currt,tFinal)
% check last N steps of time step vector: if those time steps are too small,
% we expect this to continue so that the analysis will not terminate

% init flag
abortAnalysis = false;

% remaining time
remTime = tFinal - currt;
% number of previous steps considered for criterion
N = 10;
% total steps until now
k = length(tVec);

% criterion: if sum of last N steps is smaller than a certain fraction of
%            the remaining time, abort the analysis
lastNsteps = sum(tVec(end-min(N,k)+1:end));
if remTime / lastNsteps > 1e9
    abortAnalysis = true;
    CORAwarning('CORA:contDynamics',sprintf(['The analysis is aborted because the time step size converges to 0.\n'...
        '         The reachable sets until t = ' num2str(currt) ' are returned.']));
end

end

% ------------------------------ END OF CODE ------------------------------
