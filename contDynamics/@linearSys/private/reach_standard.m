function [timeInt,timePoint,res] = reach_standard(obj,options)
% reach_standard - computes the reachable set for linear systems using the
%    standard (non-wrapping-free) reachability algorithm for linear systems
%
% Syntax:
%    [timeInt,timePoint,res] = reach_standard(obj,options)
%
% Inputs:
%    obj - linearSys object
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
%    [1] A. Girard, "Reachability of uncertain linear systems using 
%       zonotopes" in Hybrid Systems: Computation and Control, 
%       ser. LNCS 3414. Springer, 2005, pp. 291--305.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-June-2019 (from @contDynamics > reach.m)
% Last update:   19-November-2022 (MW, modularize specification check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end

%if a trajectory should be tracked
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
else
    inputCorr = 0;
end

% log information
verboseLog(1,options.tStart,options);

%initialize reachable set computations
[Rnext, options] = initReach_Euclidean(obj, options.R0, options);
Rtrans  = options.Rtrans;
Rhom    = options.Rhom;
Rhom_tp = options.Rhom_tp;
Rpar    = options.Rpar;
Raux    = options.Raux;
eAt     = obj.taylor.eAt;

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
    
    % post: ----------

    % method implemented from Algorithm 1 in [1]
    
    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,i);
        options.Rhom_tp = Rhom_tp;
        [Rhom, Rhom_tp, ~, inputCorr] = inputInducedUpdates(obj,options);
    else
        Rhom = eAt*Rhom + Rtrans;
        Rhom_tp = eAt*Rhom_tp + Rtrans;
    end
    Rhom = reduce(Rhom,options.reductionTechnique,options.zonotopeOrder);
    Rhom_tp = reduce(Rhom_tp,options.reductionTechnique,options.zonotopeOrder);
    Raux = eAt*Raux;
    Rpar = Rpar + Raux;
    if ~isnumeric(Rpar)
        Rpar = reduce(Rpar,options.reductionTechnique,options.zonotopeOrder);
    end
    
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
        if ~res
            % compute output set of last set
            timePoint.set{i+1} = outputSet(obj,options,Rnext.tp);
            return
        end
    end
    
    % log information
    verboseLog(i,tVec(i),options);
    
end

% compute output set of last set
timePoint.set{end} = outputSet(obj,options,Rstart);

% specification fulfilled
res = true;

% ------------------------------ END OF CODE ------------------------------
