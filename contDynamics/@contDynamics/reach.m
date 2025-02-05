function [R,res,options] = reach(sys,params,options,varargin)
% reach - computes the reachable continuous set for the entire time horizon
%         of a continuous system
%
% Syntax:
%    R = reach(sys,params,options)
%    [R,res] = reach(sys,params,options,spec)
%    [R,res,options] = reach(sys,params,options,spec)
%
% Inputs:
%    sys - contDynamics object
%    params - parameter defining the reachability problem
%    options - options for the computation of reachable sets
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the computed reachable set
%    res  - true if specifications are satisfied, otherwise false
%    options - options for the computation of reachable sets
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       08-August-2016
% Last update:   22-September-2016
%                28-July-2017 (EL, include algebraic reachable set)
%                20-March-2018 (NK, output sets as additional output)
%                19-May-2020 (MW, error handling for exploding sets)
%                19-November-2022 (MW, include output set computation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
spec = setDefaultValues({[]},varargin);

% options preprocessing
[params,options] = validateOptions(sys,params,options);

% handling of specifications
specLogic = [];
if ~isempty(spec)
    [spec,specLogic] = splitLogic(spec);
end

% reach called from
% - linear class (linProbSys, linParamSys) or
% - nonlinear class (nonlinearSys, nonlinDASys, nonlinParamSys)
syslin = isa(sys,'linProbSys') || isa(sys,'linParamSys');

% compute symbolic derivatives
if ~syslin
    % paramInt required for derivatives of nonlinParamSys...
    if isfield(params,'paramInt')
        options.paramInt = params.paramInt;
    end
    derivatives(sys,options);

    if contains(options.alg,'adaptive')
        % nonlinear adaptive algorithm
        [timeInt,timePoint,res,~,options] = reach_adaptive(sys,params,options);
        R = reachSet.initReachSet(timePoint,timeInt);
        return;
    end
end

% obtain factors for initial state and input solution time step
r = options.timeStep;
for i = 1:(options.taylorTerms+1)  
    options.factor(i) = r^(i)/factorial(i);    
end

% if a trajectory should be tracked
if isfield(params,'uTransVec')
    params.uTrans = params.uTransVec(:,1);
end
options.t = params.tStart;

% time period
tVec = params.tStart:options.timeStep:params.tFinal;
steps = length(tVec)-1;

% initialize cell-arrays that store the reachable set
timeInt.set = cell(steps,1);
timeInt.time = cell(steps,1);
timePoint.set = cell(steps+1,1);
timePoint.time = cell(steps+1,1);
if isa(sys,'nonlinDASys')
    timeInt.algebraic = cell(steps,1);
end

% first timePoint set is initial set
timePoint.time{1} = params.tStart;
if syslin
    timePoint.set{1} = params.R0;
else
    if isa(sys,'nonlinDASys')
        R_y = zonotope(consistentInitialState(sys,center(params.R0),...
            params.y0guess,params.uTrans));
        timePoint.set{1}{1}.set = outputSet(sys,params.R0,R_y,params,options);
    else
        timePoint.set{1}{1}.set = outputSet(sys,params.R0,params,options);
    end
    timePoint.set{1}{1}.prev = 1;
    timePoint.set{1}{1}.parent = 1;
end

% log information
verboseLog(options.verbose,1,options.t,params.tStart,params.tFinal);

% initialize reachable set computations
try
    [Rnext, options] = initReach(sys,params.R0,params,options);
catch ME
    % if error from set explosion, return corresponding information
    priv_reportReachError(ME,params.tStart,1);
    R = reachSet.initReachSet(timePoint,timeInt);
    res = false;
    return
end


% loop over all reachability steps
for i = 2:steps
    
    % save reachable set in cell structure
    if ~isa(sys,'nonlinDASys')
        timeInt.set{i-1} = outputSet(sys,Rnext.ti,params,options);
        timePoint.set{i} = outputSet(sys,Rnext.tp,params,options);
    else
        timeInt.set{i-1} = outputSet(sys,Rnext.ti,Rnext.y,params,options);
        timePoint.set{i} = outputSet(sys,Rnext.tp,Rnext.y,params,options);
    end
    timeInt.time{i-1} = interval(tVec(i-1),tVec(i));
    timePoint.time{i} = tVec(i);
    
    if isa(sys,'nonlinDASys')
        timeInt.algebraic{i-1} = Rnext.y;
    end

    % notify user if splitting has occurred
    if length(timePoint.set{i}) > length(timePoint.set{i-1})
        disp(['split! ...number of parallel sets: ',num2str(length(timePoint.set{i}))]);
    end
    
    % check specification
    if ~isempty(spec)
       if ~check(spec,Rnext.ti,timeInt.time{i-1})
           res = false;
           R = reachSet.initReachSet(timePoint,timeInt);
           return;
       end
    end

    % increment time
    options.t = tVec(i);
    % log information
    verboseLog(options.verbose,i,options.t,params.tStart,params.tFinal);

    % if a trajectory should be tracked
    if isfield(params,'uTransVec')
        params.uTrans = params.uTransVec(:,i);
    end

    % compute next reachable set
    try
        [Rnext,options] = post(sys,Rnext,params,options);
    catch ME
        % if error from set explosion, return corresponding information
        R = reachSet.initReachSet(timePoint,timeInt);
        priv_reportReachError(ME,options.t,i);
        return
    end
end

% compute output set
if ~isa(sys,'nonlinDASys')
    timeInt.set{end} = outputSet(sys,Rnext.ti,params,options);
    timePoint.set{end} = outputSet(sys,Rnext.tp,params,options);
else
    timeInt.set{end} = outputSet(sys,Rnext.ti,Rnext.y,params,options);
    timePoint.set{end} = outputSet(sys,Rnext.tp,Rnext.y,params,options);
end
timeInt.time{end} = interval(tVec(end-1),tVec(end));
timePoint.time{end} = tVec(end);

if isfield(Rnext,'y')
    timeInt.algebraic{end} = Rnext.y;
end

% check specification
if ~isempty(spec)
   if ~check(spec,timeInt.set{end},timeInt.time{end})
       res = false;
   end
end

% construct reachset object
R = reachSet.initReachSet(timePoint,timeInt);

% check temporal logic specifications
if res && ~isempty(specLogic)
    res = check(specLogic,R);
end

% log information
verboseLog(options.verbose,i+1,tVec(end),params.tStart,params.tFinal);

% ------------------------------ END OF CODE ------------------------------
