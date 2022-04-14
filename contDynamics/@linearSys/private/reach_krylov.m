function [Rout,Rout_tp,res] = reach_krylov(obj,options)
% reach_krylov - computes the reachable set for linear systems using
%  the krylov reachability algorithm for linear systems
%
% Syntax:  
%    [Rout,Rout_tp,res] = reach_krylov(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    Rout - reachable set of time intervals
%    Rout_tp - reachable set of time point
%    res  - boolean (only if specification given)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff, Mark Wetzlinger
% Written:       26-June-2019
% Last update:   14-Aug-2019
%                02-June-2020 (MA)
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end

% if a trajectory should be tracked
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
end

% log information
verboseLog(1,options.tStart,options);

%initialize reachable set computations
[Rnext,options] = initReach_Krylov(obj,options.R0,options);

%time period
tVec = options.tStart:options.timeStep:options.tFinal;

%create options.t
options.t = options.timeStep;

% initialize parameter for the output equation
Rout = cell(length(tVec)-1,1);
Rout_tp = cell(length(tVec)-1,1);


% loop over all reachability steps
for i=2:length(tVec)-1
    
    % calculate the set of system outputs   
%     Z = Rnext.ti.Z;
%     Rout{i-1} = zonotope(C*Z) + D * (options.uTrans + options.U) + k;
%     Rout{i-1} = reduce(Rout{i-1},options.reductionTechnique,options.zonotopeOrder);
    Rout{i-1} = Rnext.ti;
%     Z_tp = Rnext.tp.Z;
%     Rout_tp{i-1} = zonotope(C*Z_tp) + D * (options.uTrans + options.U) + k;
%     Rout_tp{i-1} = reduce(Rout_tp{i-1},options.reductionTechnique,options.zonotopeOrder);
    Rout_tp{i-1} = Rnext.tp;
    
    % safety property check
    if isfield(options,'specification')
        if ~check(options.specification,Rout{i-1})
            % violation
            Rout = Rout(1:i-1);
            Rout_tp = Rout_tp(1:i-1);
            res = 0;
            return
        end
    end
    
    % log information
    verboseLog(i,tVec(i),options);
    
    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,i);
    end
    
    %compute next reachable set
    [Rnext,options] = post_Krylov(obj,options);
    
    options.t = options.t + options.timeStep;
end


% Z = Rnext.ti.Z;
% Rout{end} = zonotope(C*Z) + D * (options.uTrans + options.U) + k;
% Rout{end} = reduce(Rout{end},options.reductionTechnique,options.zonotopeOrder);
Rout{end} = Rnext.ti;

% Z_tp = Rnext.tp.Z;
% Rout_tp{end} = zonotope(C*Z_tp) + D * (options.uTrans + options.U) + k;
% Rout_tp{end} = reduce(Rout_tp{end},options.reductionTechnique,options.zonotopeOrder);
Rout_tp{end} = Rnext.tp;

% safety property check
if isfield(options,'specification')
    if ~check(options.specification,Rout{i-1})
        % violation, no reduction in cell size of Rout, Rout_tp
        res = 0;
        return
    end
end

% log information
verboseLog(length(tVec),tVec(end),options);

% specification fulfilled
res = true;


end

%------------- END OF CODE --------------