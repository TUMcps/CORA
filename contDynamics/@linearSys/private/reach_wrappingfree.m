function [Rout,Rout_tp,res] = reach_wrappingfree(obj,options)
% reach - computes the reachable set for linear systems using
%  the standard (non-wrapping-free) reachability algorithm for linear systems
%
% Syntax:  
%    [Rout,Rout_tp,res] = reach_standard(obj,options)
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
% References:
%    [1] A. Girard, C. Le Guernic, and O. Maler,
%        "Efficient computation of reachable sets of
%        linear time-invariant systems with inputs"
%        in Hybrid Systems: Computation and Control, ser. LNCS 3927.
%        Springer, 2006, pp. 257--271.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff, Mark Wetzlinger
% Written:       26-June-2019 (from @contDynamics > reach.m)
% Last update:   14-Aug-2019
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i)= options.timeStep^(i)/factorial(i);    
end

% if a trajectory should be tracked
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
else
    inputCorr = 0;
end

% log information
verboseLog(1,options.tStart,options);

% initialize reachable set computations
[Rnext, options] = initReach_Euclidean(obj, options.R0, options);
Rtrans = options.Rtrans;
Rhom   = options.Rhom;
Rhom_tp = options.Rhom_tp;
Rpar   = options.Rpar;
Raux   = options.Raux;
eAt    = obj.taylor.eAt;

%time period
tVec = options.tStart:options.timeStep:options.tFinal;

% initialize parameter for the output equation
[C,D,k] = initOutputEquation(obj,options);
Rout = cell(length(tVec)-1,1);
Rout_tp = cell(length(tVec)-1,1);


% loop over all reachability steps
for i = 2:length(tVec)-1
    
    % compute output set
    [Rout{i-1},Rout_tp{i-1}] = outputSet(C,D,k,Rnext,options);
    
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
    
    % post: (wrapping free) ----------
    
    % log information
    verboseLog(i,tVec(i),options);
    
    % method implemented from Algorithm 2 in [1]
    
    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,i);
        options.Rhom_tp = Rhom_tp;
        [Rhom,Rhom_tp,Rtrans,inputCorr] = inputInducedUpdates(obj,options);
    else
        Rhom = eAt*Rhom + center(Rtrans);
        Rhom_tp = eAt*Rhom_tp + center(Rtrans);
    end
    Raux = eAt*Raux;
    Rpar = Rpar + interval(Raux) + interval(Rtrans) + (-center(Rtrans));
    
    %write results to reachable set struct Rnext
    if isa(Rhom,'mptPolytope')
        Rnext.ti = Rhom + mptPolytope(Rpar) + mptPolytope(inputCorr);
        Rnext.tp = Rhom_tp + mptPolytope(Rpar);
    else
        Rnext.ti = Rhom + zonotope(Rpar) + inputCorr;
        Rnext.tp = Rhom_tp + zonotope(Rpar);
    end

end

[Rout{end},Rout_tp{end}] = outputSet(C,D,k,Rnext,options);

if isfield(options,'specification')
    if ~check(options.specification,Rout{end})
        % violation, but no reduction in cell size of Rout, Rout_tp
        res = false;
        return
    end
end

% log information
verboseLog(length(tVec),tVec(end),options);

% specification fulfilled
res = true;

end


%------------- END OF CODE --------------