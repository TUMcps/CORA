function [R,res] = reach(linsysDT,params,options,varargin)
% reach - computes the reachable set for linear discrete-time systems
%
% Syntax:
%    R = reach(linsysDT,params,options)
%    [R,res] = reach(linsysDT,params,options,spec)
%
% Inputs:
%    linsysDT - continuous system object
%    params - model parameters
%    options - options for the computation of reachable sets
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%    res - true if specifications are satisfied, false otherwise
%
% References:
%    [1] A.A. Kurzhanskiy, P. Varaiya. "Reach set computation and control synthesis for 
%        discrete-time dynamical systems with disturbances", Automatica
%        47 (7), pp. 1414-1426, 2011.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       26-June-2019
% Last update:   08-October-2019
%                23-April-2020 (MW, restructure params/options)
%                07-December-2020 (MW, fix wrong indexing)
%                16-November-2021 (MW, add disturbance W)
%                21-December-2022 (MA, backward reachable set added)
%                09-February-2023 (LL, fix wrong uTrans indexing)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% safety property check
res = true;
spec = setDefaultValues({[]},varargin);

% options preprocessing
[params,options] = validateOptions(linsysDT,params,options);

%if a trajectory should be tracked
if isfield(params,'uTransVec')
    params.uTrans = params.uTransVec(:,1);
end

% pre-compute mappings by disturbance and noise matrices (possible since W
% and V are constant over the time horizon)
params.W = linsysDT.E * params.W;
params.V = linsysDT.F * params.V;

% compute inverse of system matrix for backward reachable sets
if strcmp(options.linAlg,'backward_maxmin') || ...
   strcmp(options.linAlg,'backward_maxmin_coarse') || ...
   strcmp(options.linAlg,'backward_maxmin_RKweighted') || ...
   strcmp(options.linAlg,'backward_minmax')
    Ainv = inv(linsysDT.A);
end
    
% time period and number of steps
tVec = params.tStart:linsysDT.dt:params.tFinal;
steps = length(tVec)-1;
options.i = 0;

% initialize parameters for the output equation
Rout = cell(steps+1,1);

% initialize input set
if iscell(linsysDT.B)
    Uadd = linsysDT.B{1}*(params.U + params.uTrans);
else
    Uadd = linsysDT.B*(params.U + params.uTrans);
end

% initialize reachable set
Rnext.tp = params.R0;

% loop over all reachability steps
for i = 1:steps
    
    % step counter
    options.i = i;    

    % compute output set at beginning of current step
    Rout{i} = outputSet(linsysDT,Rnext.tp,params,options);
    
    % write results to reachable set struct Rnext
    if strcmp(options.linAlg,'standard')
        % standard forward reachable set
        if iscell(linsysDT.A)
            Rnext.tp = linsysDT.A{i}*Rnext.tp + Uadd + linsysDT.c + params.W;
        else
            Rnext.tp = linsysDT.A*Rnext.tp + Uadd + linsysDT.c + params.W;
        end
        Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);
    elseif strcmp(options.linAlg,'backward_maxmin')
        % under-approximative backward reachable set according to 
        % Theorem 2.4 of [1]
        Rnext.tp = Ainv*(minkDiff(Rnext.tp + (-1*Uadd) - linsysDT.c,params.W,'outer'));
        Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);
    elseif strcmp(options.linAlg,'backward_maxmin_coarse')
        % under-approximative backward reachable set according to 
        % Theorem 2.4 of [1]
        Rnext.tp = Ainv*(minkDiff(Rnext.tp + (-1*Uadd) - linsysDT.c,params.W,'outer:coarse'));
        Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);
    elseif strcmp(options.linAlg,'backward_minmax')
        % under-approximative backward reachable set according to 
        % Theorem 2.4 of [1]; formula in paper not 100% correct
        Rnext.tp = Ainv*(minkDiff(Rnext.tp,params.W,'inner') + (-1*Uadd) - linsysDT.c); % <-- needs to be fixed! Rnext.tp has to be replaced
        Rnext.tp = reduceUnderApprox(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);
    elseif strcmp(options.linAlg,'backward_minmax_RKweighted')
        % under-approximative backward reachable set according to 
        % Theorem 2.4 of [1]; formula in paper not 100% correct
        Rnext.tp = Ainv*(minkDiff(Rnext.tp,params.W,'inner:RaghuramanKoeln_weighted') + (-1*Uadd) - linsysDT.c); % <-- needs to be fixed! Rnext.tp has to be replaced
        Rnext.tp = reduceUnderApprox(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);
    end

    % log information
    verboseLog(options.verbose,i,tVec(i),params.tStart,params.tFinal);

    % safety property check
    if ~isempty(spec)
        if ~check(spec,Rout{i},interval(tVec(i)))
            % violation
            timePoint.set = Rout(1:i);
            timePoint.time = num2cell(tVec(1:i)');
            R = reachSet(timePoint);
            res = false;
            return
        end
    end

    % if a trajectory should be tracked
    if isfield(params,'uTransVec') && i ~= steps
        params.uTrans = params.uTransVec(:,i+1);
        % update input set
        if iscell(linsysDT.B)
            Uadd = linsysDT.B{i}*(params.U + params.uTrans);
        else
            Uadd = linsysDT.B*(params.U + params.uTrans);
        end
    end
end

% compute last output set
Rout{end} = outputSet(linsysDT,Rnext.tp,params,options);

% construct reachable set object
timePoint.set = Rout;
timePoint.time = num2cell(tVec');
R = reachSet(timePoint);

% log information
verboseLog(options.verbose,length(tVec),tVec(end),params.tStart,params.tFinal);

% ------------------------------ END OF CODE ------------------------------
