function [Rout,Rout_tp,res] = reach_fromStart(obj,options)
% reach - computes the reachable set for linear systems using the
%  propagation of the homogeneous solution from the start
%
% Syntax:  
%    [Rout,Rout_tp,res] = reach_fromStart(obj,options)
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

% Author:        Mark Wetzlinger
% Written:       26-June-2019
% Last update:   14-Aug-2019
%                16-February-2021 (MW, correct implementation of uTransVec)
% Last revision: ---

%------------- BEGIN CODE --------------

% obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    % compute initial state factor
    options.factor(i)= options.timeStep^(i)/factorial(i);    
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


% time period
tVec = options.tStart:options.timeStep:options.tFinal;

% initialize parameter for the output equation
[C,D,k] = initOutputEquation(obj,options);
Rout = cell(length(tVec)-1,1);
Rout_tp = cell(length(tVec)-1,1);


for i=2:length(tVec)-1
    
    % calculate output
    [Rout{i-1},Rout_tp{i-1}] = outputSet(C,D,k,Rnext,options);
    
    % safety property check (only time interval)
    if isfield(options,'specification')
        if ~check(options.specification,Rout{i-1})
            % violation
            Rout = Rout(1:i-1);
            Rout_tp = Rout_tp(1:i-1);
            res = false;
            return
        end
    end
    
    
    % post --------------
    
    % log information
    verboseLog(i,tVec(i),options);
    
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