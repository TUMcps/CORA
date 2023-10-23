function [timeInt,timePoint,res] = reach_krylov(obj,options)
% reach_krylov - computes the reachable set for linear systems using
% 	 the krylov reachability algorithm for linear systems [1]
%
% Syntax:
%    [timeInt,timePoint,res] = reach_krylov(obj,options)
%
% Inputs:
%    obj - linearSys object
%    options - options for the computation of reachable sets
%
% Outputs:
%    timeInt - array of time-interval output sets
%    timePoint - array of time-point output sets
%    res - true/false (only if specification given)
%
% Example: 
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
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
end

% log information
verboseLog(1,options.tStart,options);

%initialize reachable set computations
[Rnext,options] = initReach_Krylov(obj,options.R0,options);

%time period
tVec = options.tStart:options.timeStep:options.tFinal;
steps = length(tVec)-1;

%create options.t
options.t = options.timeStep;

% initialize arguments for the output equation
timeInt.set = cell(steps,1);
timeInt.time = cell(steps,1);
timePoint.set = cell(steps+1,1);
timePoint.time = num2cell(tVec');


% loop over all reachability steps
for i=2:steps
    
    % calculate the output set
    timeInt.set{i-1} = Rnext.ti;
    timeInt.time{i-1} = interval(tVec(i-1),tVec(i));
    timePoint.set{i-1} = Rnext.tp;
    
    % safety property check
    if isfield(options,'specification')
        [res,timeInt,timePoint] = checkSpecification(...
            options.specification,[],timeInt,timePoint,i-1);
        if ~res; return; end
    end
    
    % log information
    verboseLog(i,tVec(i),options);
    
    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,i);
    end
    
    %compute next reachable set
    [Rnext,options] = post_Krylov(obj,options);
    
    % increment time
    options.t = options.t + options.timeStep;
end

% last set
timeInt.set{end} = Rnext.ti;
timeInt.time{end} = interval(tVec(end-1),tVec(end));
timePoint.set{end} = Rnext.tp;

% safety property check
if isfield(options,'specification')
    [res,timeInt,timePoint] = checkSpecification(...
        options.specification,[],timeInt,timePoint,steps);
    if ~res; return; end
end

% log information
verboseLog(length(tVec),tVec(end),options);

% specification fulfilled
res = true;

% ------------------------------ END OF CODE ------------------------------
