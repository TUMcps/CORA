function [R,res] = reach(obj,params,options,varargin)
% reach - computes the reachable set for linear discrete-time systems
%
% Syntax:  
%    R = reach(obj,params,options)
%    [R,res] = reach(obj,params,options,spec)
%
% Inputs:
%    obj - continuous system object
%    params - model parameters
%    options - options for the computation of reachable sets
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%    res - true if specifications are satisfied, false otherwise
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Mark Wetzlinger
% Written:       26-June-2019
% Last update:   08-Oct-2019
%                23-April-2020 (MW, restructure params/options)
%                07-December-2020 (MW, fix wrong indexing)
%                16-November-2021 (MW, add disturbance W)
%                09-February-2023 (LL, fix wrong uTrans indexing)
% Last revision: ---


%------------- BEGIN CODE --------------

% safety property check
res = true;

% options preprocessing
options = validateOptions(obj,mfilename,params,options);

% specification
spec = [];
if nargin >= 4
   spec = varargin{1}; 
end

%if a trajectory should be tracked
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
end

% time period and number of steps
tVec = options.tStart:obj.dt:options.tFinal;
steps = length(tVec)-1;
options.i = 0;

% initialize parameters for the output equation
Rout = cell(steps+1,1);

% initialize input set
Uadd = obj.B*(options.U + options.uTrans);

% initialize reachable set
Rnext.tp = options.R0;

% loop over all reachability steps
for i = 1:steps
    
    % step counter
    options.i = i;
    
    % compute output set at beginning of current step
    Rout{i} = outputSet(obj,options,Rnext.tp);
    
    % write results to reachable set struct Rnext
    Rnext.tp = reduce(obj.A*Rnext.tp + Uadd + obj.c + options.W,...
        options.reductionTechnique,options.zonotopeOrder);

    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,i+1);
        % update input set
        Uadd = obj.B*(options.U + options.uTrans);
    end

    % log information
    verboseLog(i,tVec(i),options);

    % safety property check
    if ~isempty(spec)
        if ~check(spec,Rout{i},interval(tVec(i)))
            % violation
            timePoint.set = Rout(1:i);
            timePoint.time = num2cell(tVec(1:i)');
            R = reachSet(timePoint);
            res = false;
            return
        end
    end
end

% compute last output set
Rout{end} = outputSet(obj,options,Rnext.tp);


% construct reachable set object
timePoint.set = Rout;
timePoint.time = num2cell(tVec');
R = reachSet(timePoint);

% log information
verboseLog(length(tVec),tVec(end),options);

%------------- END OF CODE --------------