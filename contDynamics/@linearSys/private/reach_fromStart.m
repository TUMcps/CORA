function [timeInt,timePoint,res] = reach_fromStart(linsys,params,options)
% reach_fromStart - computes the reachable set for linear systems using the
%    propagation of the homogeneous and particular solution from the start
%    note that this algorithm does not support a varying input vector
%
% Syntax:
%    [timeInt,timePoint,res] = reach_fromStart(linsys,params,options)
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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       26-June-2019
% Last update:   14-August-2019
%                16-February-2021 (MW, correct implementation of uTransVec)
%                19-November-2022 (MW, modularize specification check)
%                03-April-2024 (MW, full rewrite)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% put system into canonical form
if isfield(params,'uTransVec')
    throw(CORAerror('CORA:notSupported'));
else
    [linsys,U,u,V,v] = canonicalForm(linsys,params.U,params.uTrans,...
        params.W,params.V,zeros(linsys.nrOfNoises,1));
end

% time period, number of steps, step counter
tVec = params.tStart:options.timeStep:params.tFinal;
steps = length(tVec) - 1;

% initialize output variables for reachable sets and output sets
timeInt.set = cell(steps,1);
timeInt.time = cell(steps,1);
timePoint.set = cell(steps+1,1);
timePoint.time = num2cell(tVec');

% log information
verboseLog(options.verbose,1,params.tStart,params.tStart,params.tFinal);

% compute reachable sets for first step
[Rtp,Rti,Htp_0,Hti_0,PU_0,Pu_0] = oneStep(linsys,...
    params.R0,U,u(:,1),options.timeStep,options.taylorTerms);

% read out propagation matrix and base particular solution
eAdt = getTaylor(linsys,'eAdt',struct('timeStep',options.timeStep));

PU = PU_0;
P = eye(linsys.nrOfStates);
Q = eAdt;

% compute output sets
timePoint.set{1} = outputSet_canonicalForm(linsys,params.R0,V,v,1);
timePoint.set{2} = outputSet_canonicalForm(linsys,Rtp,V,v,2);
timeInt.set{1} = outputSet_canonicalForm(linsys,Rti,V,v,1);
timeInt.time{1} = interval(tVec(1),tVec(2));

% safety property check
if isfield(options,'specification')
    [res,timeInt,timePoint] = checkSpecification(...
        options.specification,Rti,timeInt,timePoint,1);
    if ~res; return; end
end


for k = 2:steps
    
    % propagate affine solution (no reduction since representation size
    % does not increase over time)
    Hti = Q * Hti_0;
    Htp = Q * Htp_0;
    % propagate particular solution (constant input solution is one step of
    % the propagation matrix behind because it is already included in the
    % affine solution above)
    if isa(PU,'contSet')
        PU = reduce(PU + Q * PU_0 + P * Pu_0, ...
            options.reductionTechnique, options.zonotopeOrder);
    end
    
    % R([t_k, t_k + Delta t_k]) = H([t_k, t_k + Delta t_k]) + P([0, t_k])
    Rti = Hti + PU;
    Rtp = Htp + PU;

    % compute output set
    timePoint.set{k+1} = outputSet_canonicalForm(linsys,Rtp,V,v,k+1);
    timeInt.set{k} = outputSet_canonicalForm(linsys,Rti,V,v,k);
    timeInt.time{k} = interval(tVec(k),tVec(k+1));

    % propagate matrix exponentials
    P = Q;
    Q = Q * eAdt;
    
    % safety property check
    if isfield(options,'specification')
        [res,timeInt,timePoint] = checkSpecification(...
            options.specification,Rti,timeInt,timePoint,k);
        if ~res; return; end
    end
    
    % log information
    verboseLog(options.verbose,k,tVec(k),params.tStart,params.tFinal);
    
end

% specification fulfilled
res = true;

% ------------------------------ END OF CODE ------------------------------
