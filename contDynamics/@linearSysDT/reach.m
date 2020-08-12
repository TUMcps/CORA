function [R,res] = reach(obj,params,options,varargin)
% reach - computes the reachable set for linear system
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
%    res  - 1 if specifications are satisfied, 0 if not
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Mark Wetzlinger
% Written:       26-June-2019
% Last update:   08-Oct-2019
%                23-April-2020 (restructure params/options)
% Last revision: ---


%------------- BEGIN CODE --------------

    res = 1;

    % options preprocessing
    options = params2options(params,options);
    options = checkOptionsReach(obj,options);
    
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
            if ~check(options.specification,Rout{i-1})
                % violation
                timePoint.set = Rout(1:i-1);
                timePoint.time = num2cell(tVec(2:i)');
                R = reachSet(timePoint);
                res = false;
                return
            end
        end

        % if a trajectory should be tracked
        if isfield(options,'uTransVec')
            options.uTrans = options.uTransVec(:,i);
            % update input set
            Uadd = obj.B*(options.U + options.uTrans);
        end

        % write results to reachable set struct Rnext
        if isempty(obj.c)
            Rnext.tp = obj.A*Rnext.tp + Uadd;
        else
            Rnext.tp = obj.A*Rnext.tp + Uadd + obj.c;
        end
    end

    Rout{end} = outputSet(C,D,k,Rnext,options);

    % safety property check
    res = true;

    if ~isempty(spec)
        if ~check(options.specification,Rout{end})
            % violation, but no reduction in cell size of Rout
            res = false;
        end
    end

    % construct reachable set object
    timePoint.set = Rout;
    timePoint.time = num2cell(tVec(2:end)');

    R = reachSet(timePoint);

end

%------------- END OF CODE --------------