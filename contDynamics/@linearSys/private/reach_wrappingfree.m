function [timeInt,timePoint,res] = reach_wrappingfree(obj,options)
% reach_wrappingfree - computes the reachable set for linear systems using
%    the wrapping-free reachability algorithm for linear systems [1]
%
% Syntax:
%    [timeInt,timePoint,res] = reach_wrappingfree(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    timeInt - array of time-interval reachable / output sets
%    timePoint - array of time-point reachable / output sets
%    res - true/false whether specification satisfied
%
% Example:
%    -
%
% References:
%    [1] A. Girard, C. Le Guernic, and O. Maler, "Efficient computation of
%        reachable sets of linear time-invariant systems with inputs"
%        in Hybrid Systems: Computation and Control, ser. LNCS 3927.
%        Springer, 2006, pp. 257--271.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-June-2019 (from @contDynamics > reach.m)
% Last update:   14-August-2019
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
    tracking = true;
else
    inputCorr = 0;
    tracking = false;
end

% log information
verboseLog(1,options.tStart,options);

% initialize reachable set computations
[Rnext, options] = initReach_Euclidean(obj, options.R0, options);
Rtrans = options.Rtrans;
Rhom   = options.Rhom;
Rhom_tp = options.Rhom_tp;
Rpar   = options.Rpar;
Raux   = options.Raux;
eAt    = obj.taylor.eAt;

% time period and number of steps
tVec = options.tStart:options.timeStep:options.tFinal;
steps = length(tVec) - 1;
options.i = 1;

% initialize output variables for reachable sets and output sets
timeInt.set = cell(steps,1);
timeInt.time = cell(steps,1);
timePoint.set = cell(steps+1,1);
timePoint.time = num2cell(tVec');

% compute output set of first step
timePoint.set{1} = outputSet(obj,options,options.R0);
timeInt.set{1} = outputSet(obj,options,Rnext.ti);
timeInt.time{1} = interval(tVec(1),tVec(2));
Rstart = Rnext.tp;

% safety property check
if isfield(options,'specification')
    [res,timeInt,timePoint] = checkSpecification(...
        options.specification,Rnext.ti,timeInt,timePoint,1);
    if ~res; return; end
end


% loop over all reachability steps
for i = 2:steps
    
    options.i = i;
    
    % post: (wrapping free) ----------
    
    % method implemented from Algorithm 2 in [1]
    
    % if a trajectory should be tracked
    if tracking
        options.uTrans = options.uTransVec(:,i);
        options.Rhom_tp = Rhom_tp;
        [Rhom,Rhom_tp,Rtrans,inputCorr] = inputInducedUpdates(obj,options);
    else
        Rhom = eAt*Rhom + center(Rtrans);
        Rhom_tp = eAt*Rhom_tp + center(Rtrans);
    end
    Raux = eAt*Raux;
    Rpar = Rpar + interval(Raux) + interval(Rtrans) + (-center(Rtrans));
    
    %write results to reachable set struct Rnext
    if isa(Rhom,'polytope')
        Rnext.ti = Rhom + polytope(Rpar) + polytope(inputCorr);
        Rnext.tp = Rhom_tp + polytope(Rpar);
    else
        Rnext.ti = Rhom + zonotope(Rpar) + inputCorr;
        Rnext.tp = Rhom_tp + zonotope(Rpar);
    end
    
    % ----------

    % compute output set
    timePoint.set{i} = outputSet(obj,options,Rstart);
    timeInt.set{i} = outputSet(obj,options,Rnext.ti);
    timeInt.time{i} = interval(tVec(i),tVec(i+1));
    
    % save reachable set
    Rstart = Rnext.tp;
    
    % safety property check
    if isfield(options,'specification')
        [res,timeInt,timePoint] = checkSpecification(...
            options.specification,Rnext.ti,timeInt,timePoint,i);
        if ~res; return; end
    end
    
    % log information
    verboseLog(i,tVec(i),options);
    
end

% compute output set of last set
timePoint.set{end} = outputSet(obj,options,Rstart);

% specification fulfilled
res = true;

% ------------------------------ END OF CODE ------------------------------
