function [timeInt,timePoint,res] = priv_reach_krylov(linsys,params,options)
% priv_reach_krylov - computes the reachable set for linear systems using
% 	 the krylov reachability algorithm for linear systems [1]
%
% Syntax:
%    [timeInt,timePoint,res] = priv_reach_krylov(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    timeInt - array of time-interval output sets
%    timePoint - array of time-point output sets
%    res - true/false (only if specification given)
%
% Example:
%    -
%
% References:
%    [1] M. Althoff. "Reachability analysis of large linear systems with
%        uncertain inputs in the Krylov subspace", IEEE Transactions on
%        Automatic Control 65 (2), pp. 477-492, 2020.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-June-2019
% Last update:   14-August-2019
%                02-June-2020 (MA)
%                19-November-2022 (MW, modularize specification check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end

% if a trajectory should be tracked
if isfield(params,'uTransVec')
    params.uTrans = params.uTransVec(:,1);
end

% log information
verboseLog(options.verbose,1,params.tStart,params.tStart,params.tFinal);

%initialize reachable set computations
[Rnext,options] = priv_initReach_Krylov(linsys,params.R0,options);

%time period
tVec = params.tStart:options.timeStep:params.tFinal;
steps = length(tVec)-1;

%create options.t
options.t = options.timeStep;

% initialize arguments for the output equation
timeInt.set = cell(steps,1);
timeInt.time = cell(steps,1);
timePoint.set = cell(steps+1,1);
timePoint.time = num2cell(tVec');


% loop over all reachability steps
for k=2:steps
    
    % calculate the output set
    timeInt.set{k-1} = Rnext.ti;
    timeInt.time{k-1} = interval(tVec(k-1),tVec(k));
    timePoint.set{k-1} = Rnext.tp;
    
    % safety property check
    if isfield(options,'specification')
        [res,timeInt,timePoint] = priv_checkSpecification(...
            options.specification,[],timeInt,timePoint,k-1);
        if ~res; return; end
    end
    
    % log information
    verboseLog(options.verbose,k,tVec(k),params.tStart,params.tFinal);
    
    % if a trajectory should be tracked
    if isfield(params,'uTransVec')
        params.uTrans = params.uTransVec(:,k);
    end
    
    %compute next reachable set
    [Rnext,options] = priv_post_Krylov(linsys,params,options);
    
    % increment time
    options.t = options.t + options.timeStep;
end

% last set
timeInt.set{end} = Rnext.ti;
timeInt.time{end} = interval(tVec(end-1),tVec(end));
timePoint.set{end} = Rnext.tp;

% safety property check
if isfield(options,'specification')
    [res,timeInt,timePoint] = priv_checkSpecification(...
        options.specification,[],timeInt,timePoint,steps);
    if ~res; return; end
end

% log information
verboseLog(options.verbose,length(tVec),tVec(end),params.tStart,params.tFinal);

% specification fulfilled
res = true;

% ------------------------------ END OF CODE ------------------------------
