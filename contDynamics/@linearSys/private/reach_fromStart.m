function [timeInt,timePoint,res] = reach_fromStart(obj,options)
% reach_fromStart - computes the reachable set for linear systems using the
%    propagation of the homogeneous solution from the start
%
% Syntax:
%    [timeInt,timePoint,res] = reach_fromStart(obj,options)
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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       26-June-2019
% Last update:   14-August-2019
%                16-February-2021 (MW, correct implementation of uTransVec)
%                19-November-2022 (MW, modularize specification check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    % compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end

% if a trajectory should be tracked
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
end

% log information
verboseLog(1,options.tStart,options);

% init step 
[Rnext, options] = initReach_Euclidean(obj, options.R0, options);
% init for loop
Rtrans = options.Rtrans;
Rinhom = options.Rpar;
Raux = options.Raux;
if ~isfield(options,'uTransVec')
    Rhom_0 = options.Rhom;
    Rhom_tp_0 = options.Rhom_tp;
    inputCorr = 0;
else
    Rhom_tp = options.Rhom_tp;
end
eADelta = obj.taylor.eAt;
P = eye(obj.dim);
Q = obj.taylor.eAt;


% time period, number of steps, step counter
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


for i = 2:steps
    
    options.i = i;
    
    % post --------------
    
    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,i);
        options.Rhom_tp = Rhom_tp;
        % recompute Rtrans/inputCorr -> also propagate Rhom/Rhom_tp
        [Rhom,Rhom_tp,Rtrans,inputCorr] = inputInducedUpdates(obj,options);
        % reduce homogeneous part
        Rhom_tp = reduce(Rhom_tp,options.reductionTechnique,options.zonotopeOrder);
        % propagate inhomogeneous part
        Rinhom = eADelta * Rinhom + Raux;
    else
        % propagate homogeneous part (from start -> take Rhom_0)
        Rhom    = Q * Rhom_0;
        Rhom_tp = Q * Rhom_tp_0;
        % no uTransVec -> Rtrans is constant
        Rinhom = Rinhom + Q * Raux + P * Rtrans;
    end
    
    % reduce (only inhomogeneous solution)
    Rinhom = reduce(Rinhom,options.reductionTechnique,options.zonotopeOrder);
    
    % R([t_k, t_k + Delta t_k]) = H([t_k, t_k + Delta t_k]) + P([0, t_k])
    Rnext.ti = Rhom + Rinhom + inputCorr;
    Rnext.tp = Rhom_tp + Rinhom;
    
    % propagate matrix exponentials
    P = Q;
    Q = Q * eADelta;
    
    % --------------

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
