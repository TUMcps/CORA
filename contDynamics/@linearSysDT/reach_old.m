function [R,res] = reach_old(obj,params,options,varargin)
% reach - computes the reachable set for linear system
%
% Syntax:  
%    [R,res] = reach_old(obj,params,options,varargin)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    Rout - reachable set of time intervals for the continuous dynamics
%    Rout_tp - reachable set of time points for the continuous dynamics
%    res  - boolean (only if property checked)
%    tVec - vector of time steps (only adaptive)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Mark Wetzlinger
% Written:       26-June-2019
% Last update:   08-Oct-2019
% Last revision: ---


%------------- BEGIN CODE --------------

% safety property check
res = true;

% new options preprocessing
options = validateOptions(obj,'reach',params,options);

% specification
spec = [];
if nargin >= 4
   spec = varargin{1}; 
end

%if a trajectory should be tracked
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
end

%time period
tVec = options.tStart:obj.dt:options.tFinal;

% initialize parameter for the output equation
[C,D,k] = initOutputEquation(obj,options);
Rout = cell(length(tVec)-1,1);

% initialize input set
Uadd = obj.B*(options.U + options.uTrans);

% initialize reachable set
Rnext.tp = options.R0;

% loop over all reachability steps
for i = 2:length(tVec)-1

    % compute output set
    Rout{i-1} = outputSet(C,D,k,Rnext,options);
    
    % safety property check
    if ~isempty(spec)
        if ~check(spec,Rout{i})
            % violation
            timePoint.set = Rout(1:i);
            timePoint.time = num2cell(tVec(2:i+1)');
            R = reachSet(timePoint);
            res = false;
            return
        end
    end
    
    
    % post
    
    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,i);
        
        % update input set
        Uadd = obj.B*(options.U + options.uTrans);
    end
    
    

    %write results to reachable set struct Rnext
    if isempty(obj.c)
        Rnext.tp = obj.A*Rnext.tp + Uadd;
    else
        Rnext.tp = obj.A*Rnext.tp + Uadd + obj.c;
    end
    Rnext.ti = []; %discrete time model

end

Rout{end} = outputSet(C,D,k,Rnext,options);

% safety property check
if ~isempty(spec)
    if ~check(spec,Rout{end})
        % violation, but no reduction in cell size of Rout, Rout_tp
        res = false;
        return
    end
end

% construct reachable set object
timePoint.set = Rout;
timePoint.time = num2cell(tVec(2:end)');
R = reachSet(timePoint);

end


%------------- END OF CODE --------------