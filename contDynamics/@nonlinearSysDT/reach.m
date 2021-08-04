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
%    res  - 1 if specifications are satisfied, 0 if not
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  29-January-2018
% Last revision:---

%------------- BEGIN CODE --------------

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
    R = cell(steps+1,1);
    R{1} = params.R0;

    % loop over all reachablity steps
    for i = 1:steps

        % if a trajectory should be tracked
        if isfield(options,'uTransVec')
            options.uTrans = options.uTransVec(:,i);
        end  

        % compute next reachable set
        R{i+1} = linReach(obj,R{i},options);

        % log information
        verboseLog(i,t(i),options);
        
        % check specification
        if ~isempty(spec)
           if ~check(spec,R{i+1})
               timePoint.set = R(2:i+1);
               timePoint.time = num2cell(t(2:i+1)');
               R = reachSet(timePoint);
               return;
           end
        end
    end

    % create reachable set object
    timePoint.set = R(2:end);
    timePoint.time = num2cell(t(2:end)');
    R = reachSet(timePoint);
    
    % log information
    verboseLog(length(t),t(end),options);

end

%------------- END OF CODE --------------