function [Rout,Rout_tp,res] = reach_fromStart_tVec(obj,params,options,tVec)
% reach - computes the reachable set for linear systems using the
%  propagation of the homogeneous solution from the start
% note: this function exists merely for debugging purposes
%
% Syntax:  
%    [Rout,Rout_tp,res] = reach_fromStart(obj,options)
%
% Inputs:
%    obj - continuous system object
%    params - model parameters
%    options - options for the computation of reachable sets
%    tVec - vector of time step sizes
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
% Last revision: ---


%------------- BEGIN CODE --------------

% options preprocessing
options = param2options(params,options);
options = checkOptionsReach(obj,options);

% time step
options.timeStep = tVec(1);
% obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    % compute initial state factor
    options.factor(i)= options.timeStep^(i)/factorial(i);    
end

% if a trajectory should be tracked
inputCorr = 0;

% init step 
[Rnext, options] = initReach_Euclidean(obj, options.R0, options);
% init for loop
eADelta = obj.taylor.eAt;
P = eye(obj.dim);
Q = obj.taylor.eAt;
lastRtrans = options.Rtrans;
Rtrans = lastRtrans;
Rhom_0 = options.Rhom;
Rhom_tp_0 = options.Rhom_tp;
Rinhom = options.Rpar;
Raux = options.Raux;

% initialize parameter for the output equation
[C,D,k] = initOutputEquation(obj,options);
Rout = cell(length(tVec)-1,1);
Rout_tp = cell(length(tVec)-1,1);


for i=2:length(tVec)
    
    % calculate output
    [Rout{i-1},Rout_tp{i-1}] = outputSet(C,D,k,Rnext,options);
    
    % safety property check (only time interval)
    if isfield(options,'specification')
        if inUnsafeSet(Rout{i-1},options.specification)
            % violation
            Rout = Rout(1:i-1);
            Rout_tp = Rout_tp(1:i-1);
            res = false;
            return
        end
    end
    
    
    % post --------------
    
    if tVec(i) ~= options.timeStep
        % read new timestep
        options.timeStep = tVec(i);
        % new factor
        % obtain factors for initial state and input solution
        for fac=1:(options.taylorTerms+1)
            % compute initial state factor
            options.factor(fac)= options.timeStep^(fac)/factorial(fac);    
        end
        [~, options] = initReach_Euclidean(obj, options.R0, options);
        % init for loop
        eADelta = obj.taylor.eAt;
        Rtrans = options.Rtrans;
        Rhom_0 = options.Rhom;
        Rhom_tp_0 = options.Rhom_tp;
        Raux = options.Raux;
    end
    
    % propagate homogeneous part (from start -> take Rhom_0)
    Rhom    = Q * Rhom_0;
    Rhom_tp = Q * Rhom_tp_0;
    
    % propagate inhomogeneous part
    Rinhom = Rinhom + Q * Raux + P * lastRtrans + inputCorr;
    Rinhom = reduce(Rinhom,options.reductionTechnique,options.zonotopeOrder);
    lastRtrans = Rtrans; % only important if uTransVec or non-const. timeStep
    
    % R([t_k, t_k + Delta t_k]) = H([t_k, t_k + Delta t_k]) + P([0, t_k])
    Rnext.ti = Rhom + Rinhom;
    Rnext.tp = Rhom_tp + Rinhom;
    
    % propagate matrix exponentials
    P = Q;
    Q = Q * eADelta;
    
end

[Rout{length(tVec)},Rout_tp{length(tVec)}] = outputSet(C,D,k,Rnext,options);

if isfield(options,'specification')
    if inUnsafeSet(Rout{end},options.specification)
        % violation, but no reduction in cell size of Rout, Rout_tp
        res = false;
        return
    end
end

res = true;


end

%------------- END OF CODE --------------