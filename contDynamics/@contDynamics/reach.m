function [R,res,options] = reach(obj,params,options,varargin)
% reach - computes the reachable continuous set for the entire time horizon
%         of a continuous system
%
% Syntax:
%    R = reach(obj,params,options)
%    [R,res] = reach(obj,params,options,spec)
%    [R,res,options] = reach(obj,params,options,spec)
%
% Inputs:
%    obj - contDynamics object
%    params - parameter defining the reachability problem
%    options - options for the computation of reachable sets
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the computed reachable set
%    res  - true if specifications are satisfied, otherwise false
%    options - options for the computation of reachable sets
%
% Example: 
%    
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

% options preprocessing
if ~isfield(options,'validated') || ~options.validated
    options = validateOptions(obj,mfilename,params,options);
else
    % internal call: skip validation (options already checked)
    options = validateOptions(obj,mfilename,params,options,true);
end

% handling of specifications
spec = []; specLogic = [];
if nargin >= 4
   spec = varargin{1}; 
   [spec,specLogic] = splitLogic(spec);
end

% reach called from
% - linear class (linProbSys, linParamSys) or
% - nonlinear class (nonlinearSys, nonlinDASys, nonlinParamSys)
syslin = isa(obj,'linProbSys') || isa(obj,'linParamSys');

% compute symbolic derivatives
if ~syslin
    derivatives(obj,options);
    if (isa(obj,'nonlinearSys') || isa(obj,'nonlinDASys')) && contains(options.alg,'adaptive')
        % nonlinear adaptive algorithm
        [timeInt,timePoint,res,~,options] = reach_adaptive(obj,params,options);
        % construct reachset object
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
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
end
options.t = options.tStart;

% time period
tVec = options.tStart:options.timeStep:options.tFinal;

% initialize cell-arrays that store the reachable set
timeInt.set = cell(length(tVec)-1,1);
timeInt.time = cell(length(tVec)-1,1);
timePoint.set = cell(length(tVec),1);
timePoint.time = cell(length(tVec),1);
if isa(obj,'nonlinDASys')
    timeInt.algebraic = cell(length(tVec)-1,1);
end

% first timePoint set is initial set
timePoint.time{1} = options.tStart;
if syslin
    timePoint.set{1} = options.R0;
else
    if isa(obj,'nonlinDASys')
        R_y = zonotope(consistentInitialState(obj,center(options.R0),...
            options.y0guess,options.uTrans));
        timePoint.set{1}{1}.set = outputSet(obj,options,options.R0,R_y);
    else
        timePoint.set{1}{1}.set = outputSet(obj,options,options.R0);
    end
    timePoint.set{1}{1}.prev = 1;
    timePoint.set{1}{1}.parent = 1;
end

% log information
verboseLog(1,options.t,options);

% initialize reachable set computations
try
    [Rnext, options] = initReach(obj,options.R0,options);
catch ME
    % if error from set explosion, return corresponding information
    R = [];
    reportReachError(ME,options.tStart,1);
    return
end


% loop over all reachability steps
for i = 2:length(tVec)-1
    
    % save reachable set in cell structure
    if ~isa(obj,'nonlinDASys')
        timeInt.set{i-1} = outputSet(obj,options,Rnext.ti);
        timePoint.set{i} = outputSet(obj,options,Rnext.tp);
    else
        timeInt.set{i-1} = outputSet(obj,options,Rnext.ti,Rnext.y);
        timePoint.set{i} = outputSet(obj,options,Rnext.tp,Rnext.y);
    end
    timeInt.time{i-1} = interval(tVec(i-1),tVec(i));
    timePoint.time{i} = tVec(i);
    
    if isa(obj,'nonlinDASys')
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
    verboseLog(i,options.t,options);

    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,i);
    end

    % compute next reachable set
    try
        [Rnext,options] = post(obj,Rnext,options);
    catch ME
        % if error from set explosion, return corresponding information
        R = reachSet.initReachSet(timePoint,timeInt);
        reportReachError(ME,options.t,i);
        return
    end
end

% compute output set
if ~isa(obj,'nonlinDASys')
    timeInt.set{end} = outputSet(obj,options,Rnext.ti);
    timePoint.set{end} = outputSet(obj,options,Rnext.tp);
else
    timeInt.set{end} = outputSet(obj,options,Rnext.ti,Rnext.y);
    timePoint.set{end} = outputSet(obj,options,Rnext.tp,Rnext.y);
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
verboseLog(i+1,tVec(end),options);

% ------------------------------ END OF CODE ------------------------------
