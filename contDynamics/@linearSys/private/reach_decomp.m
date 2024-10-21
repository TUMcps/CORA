function [timeInt,timePoint,res] = reach_decomp(linsys,params,options)
% reach_decomp - implementation of decomposed approach for reachability
%    analysis of linear systems, cf. [1]
%
% Syntax:
%    [timeInt,timePoint,res] = reach_decomp(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    timeInt - array of time-interval output sets
%    timePoint - array of time-point output sets
%    res - true/false whether specification satisfied
%
% Example:
%    -
%
% References: 
%   [1] S. Bogomolov, M. Forets, G. Frehse, A. Podelski, C. Schlling
%       "Decomposing Reach Set Computations with Low-dimensional Sets and
%        High-Dimensional Matrices"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       11-June-2019
% Last update:   14-August-2019
%                19-November-2022 (MW, modularize specification check)
%                16-October-2024 (MW, full rewrite)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%%% TODO: support sparsity somehow?

% put system into canonical form
if isfield(params,'uTransVec')
    throw(CORAerror('CORA:notSupported'));
else
    [linsys,U,u,V,v] = canonicalForm(linsys,params.U,params.uTrans,...
        params.W,params.V,zeros(linsys.nrOfNoises,1));
    V = V + v(:,1);
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
blocks = options.partition;
[Rtp,Rti,Htp_0,Hti_0,PU_0,Pu_0] = oneStep(linsys,...
    params.R0,U,u(:,1),options.timeStep,options.taylorTerms,blocks);

% read out propagation matrix and base particular solution
eAdt = getTaylor(linsys,'eAdt',struct('timeStep',options.timeStep));

PU = PU_0;
P = eye(linsys.nrOfStates);
Q = eAdt;

% compute output sets
timePoint.set{1} = outputSet_canonicalForm(linsys,params.R0,V,v,1);
timePoint.set{2} = outputSet_canonicalForm(linsys,recompose(Rtp),V,v,2);
timeInt.set{1} = outputSet_canonicalForm(linsys,recompose(Rti),V,v,1);
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
    Hti = block_mtimes(Q,Hti_0);
    Htp = block_mtimes(Q,Htp_0);
    % propagate particular solution (constant input solution is one step of
    % the propagation matrix behind because it is already included in the
    % affine solution above)
    PU_new = block_operation(@plus,block_mtimes(Q,PU_0),block_mtimes(P,Pu_0));
    PU = block_operation(@plus,PU,PU_new);
    if isa(PU,'contSet')
        PU = decompose(reduce(recompose(PU, ...
            options.reductionTechnique, options.zonotopeOrder),blocks));
    end
    
    % R([t_k, t_k + Delta t_k]) = H([t_k, t_k + Delta t_k]) + P([0, t_k])
    Rti = block_operation(@plus,Hti,PU);
    Rtp = block_operation(@plus,Htp,PU);

    % compute output sets
    timePoint.set{k+1} = outputSet_canonicalForm(linsys,recompose(Rtp),V,v,k+1);
    timeInt.set{k} = outputSet_canonicalForm(linsys,recompose(Rti),V,v,k);
    timeInt.time{k} = interval(tVec(k),tVec(k+1));

    % propagate matrix exponentials
    P = Q;
    Q = Q * eAdt;
    
    % safety property check
    if isfield(options,'specification')
        [res,timeInt,timePoint] = checkSpecification(...
            options.specification,recompose(Rti),timeInt,timePoint,k);
        if ~res; return; end
    end
    
    % log information
    verboseLog(options.verbose,k,tVec(k),params.tStart,params.tFinal);
    
end

% specification fulfilled
res = true;

% ------------------------------ END OF CODE ------------------------------
