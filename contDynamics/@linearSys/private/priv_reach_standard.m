function [timeInt,timePoint,res] = priv_reach_standard(linsys,params,options)
% priv_reach_standard - computes the reachable set for linear systems using the
%    standard (non-wrapping-free) reachability algorithm for linear systems
%
% Syntax:
%    [timeInt,timePoint,res] = priv_reach_standard(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    timeInt - array of time-interval reachable / output sets
%    timePoint - array of time-point reachable / output sets
%    res - true/false whether specification satisfied
%
% Example:
%    -
%
% References:
%    [1] A. Girard, "Reachability of uncertain linear systems using 
%       zonotopes" in Hybrid Systems: Computation and Control, 
%       ser. LNCS 3414. Springer, 2005, pp. 291--305.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-June-2019 (from @contDynamics > reach.m)
% Last update:   19-November-2022 (MW, modularize specification check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% time period and number of steps
tVec = params.tStart:options.timeStep:params.tFinal;
steps = length(tVec) - 1;

% put system into canonical form: this allows us to compute the output sets
% much more efficiently (see below)
if isfield(params,'uTransVec')
    [linsys,U,u,V,v] = canonicalForm(linsys,params.U,params.uTransVec,...
        params.W,params.V,zeros(linsys.nrOfNoises,1));
else
    [linsys,U,u,V,v] = canonicalForm(linsys,params.U,params.uTrans,...
        params.W,params.V,zeros(linsys.nrOfNoises,1));
end
% check if time-varying inputs given
isU = representsa(U,'origin');

% initialize output variables for reachable sets and output sets
timeInt.set = cell(steps,1);
timeInt.time = cell(steps,1);
timePoint.set = cell(steps+1,1);
timePoint.time = num2cell(tVec');

% log information
verboseLog(options.verbose,1,params.tStart,params.tStart,params.tFinal);

% compute reachable sets for first step
[Rtp,Rti,Htp,Hti,PU,Pu,~,C_input] = oneStep(linsys,...
    params.R0,U,u(:,1),options.timeStep,options.taylorTerms);
% read out propagation matrix and base particular solution
eAdt = getTaylor(linsys,'eAdt',struct('timeStep',options.timeStep));
if ~isU
    PU_next = PU;
end

% compute output set of start set and first time-interval solution
timePoint.set{1} = priv_outputSet_canonicalForm(linsys,params.R0,V,v,1);
timeInt.set{1} = priv_outputSet_canonicalForm(linsys,Rti,V,v,1);
timeInt.time{1} = interval(tVec(1),tVec(2));
timePoint.set{2} = priv_outputSet_canonicalForm(linsys,Rtp,V,v,2);

% safety property check
if isfield(options,'specification')
    [res,timeInt,timePoint] = priv_checkSpecification(...
        options.specification,Rti,timeInt,timePoint,1);
    if ~res; return; end
end


% loop over all reachability steps
for k = 2:steps

    % method implemented from Algorithm 1 in [1]
    
    % re-compute particular solution due to constant input if we have a
    % time-varying input trajectory, since the constant input is included
    % in our affine solution, we recompute Htp, Hti, and Pu, incl. errors
    if isfield(params,'uTransVec')
        Htp_start = Htp;
        [Htp,Pu,~,C_state,C_input] = affineSolution(...
            linsys,Htp_start,u(:,k),options.timeStep,options.taylorTerms);
        Hti = enclose(Htp_start,Htp) + C_state;
    else
        % homogeneous solution, incl. reduction
        Htp = eAdt * Htp + Pu;
        Hti = eAdt * Hti + Pu;
    end
    % reduction
    Htp = reduce(Htp,options.reductionTechnique,options.zonotopeOrder);
    Hti = reduce(Hti,options.reductionTechnique,options.zonotopeOrder);

    if ~isU
        % propagate particular solution (time-varying, centered at zero)
        PU_next = eAdt * PU_next;
        PU = reduce(PU + PU_next, options.reductionTechnique,...
            options.zonotopeOrder);
    end

    % compute reachable set
    Rti = Hti + PU + C_input;
    Rtp = Htp + PU;
    
    % compute output set
    timeInt.set{k} = priv_outputSet_canonicalForm(linsys,Rti,V,v,k);
    timeInt.time{k} = interval(tVec(k),tVec(k+1));
    
    % compute output set for start set of next step
    timePoint.set{k+1} = priv_outputSet_canonicalForm(linsys,Rtp,V,v,k+1);

    % safety property check
    if isfield(options,'specification')
        [res,timeInt,timePoint] = priv_checkSpecification(...
            options.specification,Rti,timeInt,timePoint,k);
        if ~res; return; end
    end
    
    % log information
    verboseLog(options.verbose,k,tVec(k),params.tStart,params.tFinal);
    
end

% specification fulfilled
res = true;

% ------------------------------ END OF CODE ------------------------------
