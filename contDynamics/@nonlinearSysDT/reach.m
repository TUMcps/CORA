function R = reach(obj,params,options,varargin)
% reach - computes the reachable sets of the discrete time system
%
% Syntax:
%    R = reach(obj,params,options)
%    [R,res] = reach(obj,params,options,spec)
%
% Inputs:
%    obj - nonlinearSysDT object
%    params - parameter defining the reachability problem
%    options - options for the computation of the reachable set
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%    res - true/false whether specifications are satisfied
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       21-August-2012
% Last update:   29-January-2018
%                19-November-2022 (MW, integrate output equation)
%                10-May-2023 (LL, integrate uTrans in U)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% options preprocessing
options = validateOptions(obj,mfilename,params,options);

spec = [];
if nargin >= 4
   spec = varargin{1}; 
end

% compute symbolic derivatives
derivatives(obj,options);

% initialize cell array that stores the reachable sets
t = options.tStart:obj.dt:options.tFinal;

steps = length(t)-1;
timePoint.set = cell(steps+1,1);

% add constant input
if isfield(options,'uTrans')
    options.U = options.U + options.uTrans;
    options.uTrans = 0;
end
U0 = options.U;

% add input for time 1
if isfield(options,'uTransVec')
    options.U = U0 + options.uTransVec(:,1);
end  

% compute output for time 1
timePoint.set{1} = outputSet(obj,options,params.R0);
Rnext = params.R0;

% loop over all reachablity steps
for i = 1:steps

    options.i = i;
    % compute next reachable set
    [Rnext,options] = linReach(obj,Rnext,options);

    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.U = U0 + options.uTransVec(:,i+1);
    end  
    % compute output set
    timePoint.set{i+1} = outputSet(obj,options,Rnext);

    % log information
    verboseLog(i,t(i),options);
    
    % check specification
    if ~isempty(spec)
       if ~check(spec,timePoint.set{i+1},interval(t(i+1)))
           timePoint.set = timePoint(1:i+1);
           timePoint.time = num2cell(t(1:i+1)');
           R = reachSet(timePoint);
           return;
       end
    end
end

% create reachable set object
timePoint.time = num2cell(t(1:end)');
R = reachSet(timePoint);

% log information
verboseLog(length(t),t(end),options);

% ------------------------------ END OF CODE ------------------------------
