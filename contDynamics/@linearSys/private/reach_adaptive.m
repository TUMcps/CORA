function [timeInt,timePoint,res,savedata] = reach_adaptive(linsys,params,options)
% reach_adaptive - computes an outer approximation of the reachable set
%    for linear time-invariant systems given a maximum error in terms of
%    the Hausdorff distance to the exact reachable set; all internal
%    reachability settings are set automatically
%
% Syntax:
%    [timeInt,timePoint,res,savedata] = reach_adaptive(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    timeInt - array of time-interval reachable / output sets
%    timePoint - array of time-point reachable / output sets
%    res - specifications verified (only if called from verify)
%    savedata - data used for subsequent runs (only if called from verify)
%
% Example:
%    -
%
% References:
%    [1] M. Wetzlinger et al. "Fully automated verification of linear
%        systems using inner-and outer-approximations of reachable sets",
%        TAC, 2023.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       08-July-2021
% Last update:   22-March-2022
% Last revision: 06-November-2024 (MW, full refactor)

% ------------------------------ BEGIN CODE -------------------------------

% initializations ---------------------------------------------------------

% call from verify?
if ~isfield(options,'verify')
    options.verify = false;
end
% for safety, init taylor helper class if not there already
if ~isa(linsys.taylor,'taylorLinSys')
    linsys.taylor = taylorLinSys(linsys.A);
end

% init struct which keeps sets and timeStep to avoid recomputations
savedata = aux_initSavedata(options);

% time intervals where specifications are not yet verified
[timeSpecUnsat,params.tFinal] = aux_initSpecUnsat(params,options.verify);

% init step and time (shift so that start at t=0), and spec satisfaction
k = 0; t = 0; res = true;
params.tFinal = params.tFinal - params.tStart; params.tStart = 0;

% initialize first time step by time horizon
timeStep = params.tFinal;

% initialize all variables for sets
set = aux_initSets(linsys,params,options.verify);

% init over-approximation error w.r.t the corresponding exact reachable
% sets (this is not a return value, but might still be convenient to have)
Rcont_error = [];
Rcont_tp_error = [];

% initialize time-point and time-interval output sets and
% over-approximation error w.r.t the corresponding exact output sets
timeInt.set = cell(0);
timeInt.time = cell(0);
timeInt.error = [];
timePoint.set = cell(0);
timePoint.time = cell(0);
timePoint.error = [];

% rewrite equations in canonical form
[linsys,params] = aux_canonicalForm(linsys,params);

% init all auxiliary flags
[isU,G_U,isu,isuconst] = aux_initFlags(params);

% compute factor for scaling of Hausdorff distance from state space to
% output space due to linear transformation with C, see [1, (44)]
errR2Y = norm(sqrtm(full(linsys.C * linsys.C')),2);  % same as norm(linsys.C)
% scale maximum error so that it corresponds to R (C = 1 => errR2Y = 1)
errs = linErrorBound(options.error/errR2Y, params.tFinal);
% time step adaptation using approximation functions or bisection:
% first step always bisection (more robust)
errs.useApproxFun = false;

% pre-compute how much of the total to allocate to the reduction error
[errs,savedata] = aux_errorboundreduction(errs,linsys.A,isU,G_U,options,savedata);

% compute first output solution and corresponding error
timePoint.set{1} = linsys.C * params.R0 + params.V + params.vTrans;
timePoint.error(1) = 0;
% -------------------------------------------------------------------------


% main loop until the time is up
while params.tFinal - t > 1e-9
    
    % update step counter
    k = k + 1;
    
    % log information
    verboseLog(options.verbose,k,t,params.tStart,params.tFinal);
    
    % update constant input vector (if necessary)
    [params.uTrans, params.vTrans] = aux_uTrans_vTrans(t,...
        params.uTransVec,params.vTransVec,params.tu);
    
    % flag whether full computation required (only where specifications are 
    % not yet satisfied), additionally return time until next change
    if numIntervals(timeSpecUnsat) == 0
        maxTimeStepSpec = Inf; fullcomp(k,1) = true;
    else
        [maxTimeStepSpec,fullcomp(k,1)] = timeUntilSwitch(timeSpecUnsat,t);
    end

    % determine startset for current step
    if fullcomp(k)
        set.startset = set.Hstartp;
    end
    
    % compute maximum admissible time step size for current step (either
    % until the time horizon or the delayed start, or shorter if there is a
    % time-varying input vector)
    maxTimeStep = aux_maxTimeStep(t,params.uTransVec,params.tu,...
        params.tFinal,maxTimeStepSpec);

    % initial guess for time step size of the current step:
    % same time step size as in the previous step (unless this exceeds the
    % value for the maximum admissible timeStep)
    timeStep = min([maxTimeStep,timeStep]);
    
    % initialize current step struct and iteration counter
    k_iter = 0;
    set = aux_prepareSets(set);
    
    % loop until both errors are satisfied
    while true
        
        % increase counter for time step sizes used in this step
        k_iter = k_iter + 1;
        
        % compute error bounds
        nextBounds(errs,timeStep,t,k,k_iter);
        
        % check if time step size was used before (index in savedata struct)
        timeStepIdx(k_iter,1) = aux_checkTimeStep(timeStep,savedata);
        
        % compute errors for given time step size
        [errs,set,converged] = aux_errors(linsys,k,k_iter,errs,fullcomp(k),...
            isU,isuconst,G_U,params.uTrans,set,...
            t,timeStep,timeStepIdx(k_iter,1),savedata);
        if ~converged
            % if computation of correction matrices and/or particular
            % solution has not converged, repeat with smaller time step
            % size and reset counter
            if k_iter == 1
                timeStep = timeStep * 0.1; k_iter = 0;
                set = aux_prepareSets(set);
                continue
            else
                % take sets for time step from previous iteration (must be ok)
                k_iter_chosen = find(timeStepIdx(1:k_iter) == timeStepIdx(k_iter-1),1,'first');
                break
            end
        end

        % if new time step size used, save auxiliary variables so that they
        % can be read in another step (or entire run, in verify-algorithm)
        if timeStepIdx(k_iter) == 0 || ...
                fullcomp(k) ~= savedata.fullcomp(timeStepIdx(k_iter))
            [timeStepIdx(k_iter,1),savedata] = aux_savedata(...
                savedata,k,k_iter,fullcomp(k),isU,isuconst,set,errs,timeStep);
        end

        % init/update coefficients of approximation functions
        updateCoefficientsApproxFun(errs,k,k_iter,fullcomp(k));
        
        % check non-accumulating and accumulating error against bound
        errs.bound_acc_ok{k,1}(k_iter,1) = errs.step_acc{k}(k_iter) < errs.bound_acc{k}(k_iter);
        errs.bound_nonacc_ok{k,1}(k_iter,1) = errs.seq_nonacc{k}(k_iter) < errs.bound_rem{k}(k_iter) || ~fullcomp(k);
        
        % update bisect values (lower/upper bounds also used in adaptation
        % approximation using functions)
        updateBisection(errs,k,k_iter,isU,timeStep);
        
        % check whether suitable time step size found
        if errs.bound_acc_ok{k}(k_iter,1) && (~fullcomp(k) || errs.bound_nonacc_ok{k}(k_iter,1))
            
            % proceed if timeStep is already maximum admissible time step
            if abs(timeStep - maxTimeStep) < 1e-12
                k_iter_chosen = k_iter;
                break;
            end
            
            % ensure that either the joint bound for eacc and enonacc
            % or the bound for eacc is fulfilled relatively closely
            e_acc_close = errs.step_acc{k}(k_iter) > 0.90 * errs.bound_acc{k}(k_iter);
            e_rem_close = errs.step_acc{k}(k_iter) + ...
                  errs.seq_nonacc{k}(k_iter) - errs.idv_PUtkplus1{k}(k_iter) ...
                > 0.90 * errs.bound_rem{k}(k_iter);
            if e_acc_close || (fullcomp(k) && e_rem_close)
                k_iter_chosen = k_iter;
                break;
            end
        end
        
        % predict timeStep necessary to satisfy both (non-acc/acc) bounds
        timeStep = estimateTimeStepSize(errs,t,k,k_iter,fullcomp(k),timeStep,maxTimeStep,isU);
        
        % adjust chosen time step size to avoid recomputations
        [timeStep,timeStepIdx(k_iter+1,1)] = ...
            aux_adjustTimeStep(k,timeStep,isU,maxTimeStep,savedata);
        
        % exit if time step size already proposed in a prior iteration
        timeStepIdx_ = find(timeStepIdx(1:k_iter) == timeStepIdx(k_iter+1),1,'first');
        if ~isempty(timeStepIdx_) && errs.bound_acc_ok{k}(timeStepIdx_) ...
                && (~fullcomp(k) || errs.bound_nonacc_ok{k}(timeStepIdx_))
            k_iter_chosen = timeStepIdx_;
            break;
        end
        
    end

    % compute propagation matrix
    eAdtk = readFieldForTimeStep(linsys.taylor,'eAdt',timeStep);
    if isU
        eAtkplus1 = eAdtk*readFieldForTimeStep(linsys.taylor,'eAdt',t);
        insertFieldTimeStep(linsys.taylor,'eAdt',eAtkplus1,t+timeStep);
    elseif ~fullcomp(k) && abs(timeStep - maxTimeStepSpec) < 1e-9
        % switch from partial computation to full computation in next step
        eAtkplus1 = computeField(linsys.taylor,'eAdt',struct('timeStep',t+timeStep));
    end

    % clean up: remove sets for other time step sizes
    set = aux_postCleanUpSets(set,k_iter_chosen,fullcomp(k),isU);
    
    % accumulate particular solution
    [set,errs.step_red(k,1)] = aux_reduce(linsys,t,set,errs.bound_red{k}(k_iter_chosen),isU);
    
    % accumulate accumulating error and reduction error, save non-acc error
    accumulateErrors(errs,k,k_iter_chosen);
    % clean up: values for other time step sizes no longer needed
    removeRedundantValues(errs,k,k_iter_chosen);
    
    % compute particular solution due to input vector
    if isu
        % due to current method of computation, this induces no error,
        % as the solution is just a point!
        if ~isuconst || isempty(savedata.Pu{timeStepIdx(k_iter_chosen)})
            set.Pu = particularSolution_constant(linsys,params.uTrans,timeStep,Inf);
            if isuconst
                savedata.Pu{timeStepIdx(k_iter_chosen),1} = set.Pu;
            end
        else
            % spare re-computation if time step size has been used before
            set.Pu = savedata.Pu{timeStepIdx(k_iter_chosen)};
        end
        
        % propagate constant input solution (first step: Putotal = 0)
        set.Putotal = eAdtk * set.Putotal + set.Pu;
    end
    
    % propagate reachable sets
    if ~fullcomp(k)
        % no computation of the full affine solution or reachable sets or 
        % output sets (incl. the contained errors)
        
        % save NaN value for errors since sets have not been computed
        Rcont_error(k,1) = NaN;
        Rcont_tp_error(k,1) = NaN;
        
        % compute time-point affine solution at the end of the step, which
        % is required to start the full computation in the next step
        if abs(timeStep - maxTimeStepSpec) < 1e-9
            set.Hstartp = eAtkplus1*params.R0 + set.Putotal;
        end
        
        % set Rout entry to empty and error to NaN for consistency with tVec
        timeInt.set{k,1} = [];
        timeInt.error(k,1) = NaN;
        
        % propagate auxiliary term for quick check of inner-approximations
        if options.verify
            [res,set,params] = aux_quickCheck(linsys,set,[t,t+timeStep],[],...
                errs.cum_red(k),Rcont_error(k),errR2Y,params);
        end
    
    else
        % full propagation of sets and errors
            
        % compute auxiliary solution of affine system
        set.Hstartp = eAdtk * set.startset + set.Pu;
        set.enc = enclose(set.startset, set.Hstartp);
        
        % gather over-approximation errors
        [Rcont_error(k,1),Rcont_tp_error(k,1)] = fullErrors(errs,k);

        % compute output set and contained over-approximation error
        [timeInt.set{k,1},timeInt.error(k,1),timePoint.set{k+1,1},timePoint.error(k+1,1),set] = ...
            aux_outputset(linsys,isU,set,Rcont_error(k),params,errR2Y,Rcont_tp_error(k));
        if options.verify
            % check violation of inner-approximation for quick exit
            [res,set,params] = aux_quickCheck(linsys,set,[t,t+timeStep],...
                center(timeInt.set{k}),errs.cum_red(k),timeInt.error(k),errR2Y,params);
        end
    
    end
    
    % save chosen time step to time step vector
    tVec(k,1) = timeStep;
    
    % update starting time of next step (needs to be done before checking
    % the loop condition with respect to the time horizon)
    t = t + timeStep;

    % time step adaptation using approximation functions or bisection:
    errs.useApproxFun = true;
    
    % exit if simple-to-check specifications are violated
    if ~res
        break
    end
    
end

% for debugging/investigating
% print(errs);
% plot(errs);

% log information
verboseLog(options.verbose,k,t,params.tStart,params.tFinal);

% check if all errors below respective bounds
if ~checkErrors(errs,Rcont_error,Rcont_tp_error)
    throw(CORAerror('CORA:reportToDev','Error bounds not respected.'));
end

% write time to output structs
timePoint.time = num2cell([0;cumsum(tVec)]+params.tStart);
timeInt.time = cell(length(timePoint.time)-1,1);
for k=1:length(timePoint.time)-1
    timeInt.time{k} = interval(timePoint.time{k},timePoint.time{k+1});
end

% save changes to specification in savedata
if options.verify
    savedata.safeSet = params.safeSet;
    savedata.unsafeSet = params.unsafeSet;
end

end


% Auxiliary functions -----------------------------------------------------

function [isU,G_U,isu,isuconst] = aux_initFlags(params)
% saving of operations for affine systems (u = const. over entire time
% horizon) vs. system with varying u or even uncertainty U
% -> the resulting if-else branching looks quite ugly, but still yields
% large speed-ups for high-dimensional systems
G_U = generators(params.U);
isU = ~isempty(G_U);
isu = isfield(params,'uTransVec') && any(params.uTransVec,'all');
isuconst = ~(isfield(params,'uTransVec') && size(params.uTransVec,2) > 1);

end

function set = aux_initSets(sys,params,verify)
% init set-struct, which comprises auxiliary reachable sets, generator
% matrices of (partial) particular solutions, and information about the
% fast inner-approximation check

% initialize auxiliary solutions
set.startset = params.R0;
set.Hstarti = [];
set.Hstartp = params.R0;
set.Pu = zeros(sys.nrOfStates,1);
set.Putotal = zeros(sys.nrOfStates,1);
set.PU = [];
set.PUtotal = zeros(sys.nrOfStates,1);
set.Rcont_tp_kminus1 = params.R0;

% for reduction and propagation of particular solution
set.girard_PUtotal_zero = [];
set.G_PUtotal_infty = zeros(sys.nrOfStates,1);
set.G_PUtotal_zero = double.empty(sys.nrOfStates,0);
set.Gabs_PUtotal_zero = double.empty(sys.nrOfStates,0); % only omega_rad
set.Ghat_PUtotal_zero = double.empty(sys.nrOfStates,0); % only omega_max
set.nrG_PUtotal_zero = 0;

% only fast inner-approximation check
if verify
    set.G_GfastInner = repmat({0},numel(params.safeSet),1);
    set.F_GfastInner = repmat({0},numel(params.unsafeSet),1);
    set.GfastInner_add = zeros(sys.nrOfOutputs,1);
end

end

function set = aux_prepareSets(set)
% init/remove previous auxiliary variables for computation of errors and
% output sets

% affine solution
set.boxFstartset_center = cell(1,1);
set.boxFstartset_Gbox = cell(1,1);
set.GuTrans_center = cell(1,1);
set.GuTrans_Gbox = cell(1,1);
% particular solution
set.G_PU_zero = cell(1,1);
set.G_PU_infty = cell(1,1);
set.PU_A_sum_error = cell(1,1);

end

function set = aux_postCleanUpSets(set,k_iter,fullcomp,isU)
% after the time step size has been determined, we remove all other
% computed auxiliary sets for further computations

if fullcomp
    set.boxFstartset_center = set.boxFstartset_center{k_iter};
    set.boxFstartset_Gbox = set.boxFstartset_Gbox{k_iter};
    set.GuTrans_center = set.GuTrans_center{k_iter};
    set.GuTrans_Gbox = set.GuTrans_Gbox{k_iter};
end
if isU
    set.G_PU_zero = set.G_PU_zero{k_iter};
    set.G_PU_infty = set.G_PU_infty{k_iter};
    set.PU_A_sum_error = set.PU_A_sum_error{k_iter};
end

end

function [linsys,params] = aux_canonicalForm(linsys,params)
% convert system to canonical form to limit if-else branching

% all good if system is already in canonical form
if isCanonicalForm(linsys)
    return
end

% initialize input vector if sequence given
if isfield(params,'uTransVec')
    % time-varying input vector
    uVec = params.uTransVec;
else
    % no time-varying input vector, but uTrans given
    uVec = params.uTrans;
end

% rewrite in canonical form
[linsys,U_,u_,V_,v_] = canonicalForm(linsys,...
    params.U,uVec,params.W,params.V,zeros(linsys.nrOfNoises,1));

params.U = U_;
params.uTransVec = u_;
params.V = V_;
params.vTransVec = v_;

% assign first values
params.uTrans = params.uTransVec(:,1);
params.vTrans = params.vTransVec(:,1);

end

function [errs,savedata] = aux_errorboundreduction(errs,A,isU,G_U,options,savedata)
% compute the percentage of the total error allocated to the reduction
% error

if ~isU
    % no reduction, entire error margin is available for nonacc errors
    errs.bound_red_max = 0;
    savedata.reductionerrorpercentage = errs.bound_red_max / errs.emax;
    return
elseif isfield(options,'reductionerrorpercentage')
    % define error curve manually
    errs.bound_red_max = errs.emax * options.reductionerrorpercentage;
    savedata.reductionerrorpercentage = errs.bound_red_max/errs.emax;
    return
elseif isfield(savedata,'reductionerrorpercentage')
    % reduction error defined by previous run (only verify)
    errs.bound_red_max = errs.emax * options.savedata.reductionerrorpercentage;
    return
end

% post-step: vector over time for reduction error allocation
computeErrorBoundReduction(errs,A,G_U);

% save value
savedata.reductionerrorpercentage = errs.bound_red_max/errs.emax;

end

function savedata = aux_initSavedata(options)
% initialize all fields in struct savedata

% read savedata from previous run
if isfield(options,'savedata')
    savedata = options.savedata;
else
    savedata = struct();
end

fieldsEmptyVector = {'timeStep','fullcomp','eG'};
for i=1:numel(fieldsEmptyVector)
    if ~isfield(savedata,fieldsEmptyVector{i})
        savedata.(fieldsEmptyVector{i}) = [];
    end
end

fieldsEmptyCell = {...
    'F','G','GuTrans_center','GuTrans_Gbox','Pu',... % affine dynamics
    'G_PU_zero','G_PU_infty','PU_A_sum_error' % particular solution
    };
for i=1:numel(fieldsEmptyCell)
    if ~isfield(savedata,fieldsEmptyCell{i})
        savedata.(fieldsEmptyCell{i}) = {};
    end
end

end


% specification-related functions
function [timeSpecUnsat,tFinal] = aux_initSpecUnsat(params,verify)
% initialize time intervals where all specifications are (not) satisfied
% note: tStart must not be shifted to 0 before calling this function!

tFinal = params.tFinal;
timeSpecUnsat = verifyTime();

if verify
    times = [cellfun(@(S) S.time_,params.unsafeSet,'UniformOutput',false); ...
             cellfun(@(S) S.time_,params.safeSet,'UniformOutput',false)];
    % shift by start time
    times = shift([times{:}]',-params.tStart);
    timeSpecUnsat = unify(times);
    tFinal = finalTime(timeSpecUnsat) + params.tStart;
end

end

function [res,set,params] = aux_quickCheck(linsys,set,timeInterval,...
    Rout_center,eredtotal,Rout_error,errR2Y,params)

% default value: all ok
res = true;
return  % currently disabled (messes saving of sets up)

% check for simple exit (safe set violated)
for i=1:numel(params.safeSet)
    
    % propagate generator matrix of GfastInner (for later steps)
    if params.safeSet{i}.fastInner
        C = params.safeSet{i}.set.A;
        % note that G_GfastInner contains mapping by C from i-th safe set,
        % but GfastInner_add does not (so it can be used for all safe sets)
        set.G_GfastInner{i} = set.G_GfastInner{i} + sum(abs(C * set.GfastInner_add),2);
        
        % perform actual check
        if isIntersecting(params.safeSet{i}.time_,timeInterval)
            d = params.safeSet{i}.set.b;
            Grem = sum(abs([C * linsys.C * generators(set.enc), ...
                C * linsys.C * diag(set.boxFstartset_Gbox + set.GuTrans_Gbox + set.G_PUtotal_infty)]),2);
            
            % previously: max(-d + C*c + sum(abs(C*G),2)) > Rout_error(k)
            % but now sum(abs(C*G),2) aims to avoid re-computations
            checkval = max(-d + C*Rout_center + set.G_GfastInner{i} + Grem);
            if checkval > Rout_error - errR2Y * eredtotal
                % inner-approximation will be outside of safe set
                res = false; return
            elseif checkval < 0
                % outer-approximation inside safe set -> time interval ok
                params.safeSet{i}.time_ = ...
                    setdiff(params.safeSet{i}.time_,timeInterval);
            end
            
        end
        
    end
    
end

% check for simple exit (unsafe set violated)
for i=1:numel(params.unsafeSet)
    
    % propagate generator matrix of GfastInner (for later steps)
    if params.unsafeSet{i}.fastInner
        C = params.unsafeSet{i}.set.A;
        % note that G_GfastInner contains mapping by C from i-th safe set,
        % but GfastInner_add does not (so it can be used for all safe sets)
        set.F_GfastInner{i} = set.F_GfastInner{i} + sum(abs(C * set.GfastInner_add),2);
        
        % perform actual check
        if isIntersecting(params.unsafeSet{i}.time_,timeInterval)
            d = params.unsafeSet{i}.set.b;
            Grem = sum(abs([C * linsys.C * generators(set.enc), ...
                C * linsys.C * diag(set.boxFstartset_Gbox + set.GuTrans_Gbox + set.G_PUtotal_infty)]),2);
            
            % previously: max(-d + C*c + sum(abs(C*G),2)) > Rout_error(k)
            % but now sum(abs(C*G),2) aims to avoid re-computations
            checkval = max(-d + C*Rout_center + set.G_GfastInner{i} + Grem);
            if checkval > Rout_error - errR2Y * eredtotal
                % inner-approximation will intersect unsafe set
                res = false; return
            elseif checkval < 0
                % outer-approximation does not intersect unsafe set ->
                % time interval is verified
                params.unsafeSet{i}.time_ = ...
                    setdiff(params.unsafeSet{i}.time_,timeInterval);
            end
            
        end
        
    end
    
end

end


% set propagation
function [uTrans,vTrans] = aux_uTrans_vTrans(t,uTransVec,vTransVec,tu)
% reads the value of the input vector (uTrans, vTrans) from the matrix
% storing all input vectors (uTransVec, vTransVec) for the current step
% (note: it is already ensured that the input vector changes over time)

% piecewise-constant offset in input vector for state equation
if size(uTransVec,2) == 1
    % constant offset
    uTrans = uTransVec(:,1);
else
    % take correct matrix index depending on current time: note that we
    % compute the output set for R([t_k,t_k+1]) and R(t_k+1), the latter of
    % which requires already the next entry of vTransVec
    idx = find(t >= tu,1,'last');
    uTrans = uTransVec(:,idx);
end

% piecewise-constant offset in input vector for output equation
if size(vTransVec,2) == 1
    % constant offset
    vTrans = vTransVec(:,1);
else
    % take correct matrix index depending on current time: note that we
    % compute the output set for R([t_k,t_k+1]) and R(t_k+1), the latter of
    % which requires already the next entry of vTransVec
    idx = find(t >= tu,1,'last');
    vTrans = vTransVec(:,idx+1);
end

end

function [Rout,Rout_error,Rout_tp,Rout_tp_error,set] = ...
	aux_outputset(linsys,isU,set,Rcont_error,params,errR2Y,Rcont_tp_error)
% computes the time-point and time-interval output set as well as the
% contained over-approximation error; note: homogeneous and particular
% solutions are first mapped by the output equation and only then added
% (saving memory allocation when the state dimension is much larger than
% the output dimension)

if isscalar(linsys.C) && linsys.C == 1 && representsa_(params.V,'origin',eps) && ~any(params.vTrans)
    % y = x ... consequently, errors are also equal
    if isU
        Rout = zonotope(...
            set.enc.c + set.boxFstartset_center + set.GuTrans_center, ...
            [set.enc.G, set.G_PUtotal_zero, ...
            diag(set.boxFstartset_Gbox + set.GuTrans_Gbox + set.G_PUtotal_infty)]);
    else
        Rout = zonotope(...
            set.enc.c + set.boxFstartset_center + set.GuTrans_center, ...
            [set.enc.G, diag(set.boxFstartset_Gbox + set.GuTrans_Gbox)]);
    end
    Rout_error = Rcont_error;
    
    % time-point solution
    if isU
        Rout_tp = zonotope(set.Hstartp.c, ...
            [set.Hstartp.G, set.G_PUtotal_zero, diag(set.G_PUtotal_infty)]);
    else
        Rout_tp = set.Hstartp;
    end
    Rout_tp_error = Rcont_tp_error;
    
elseif ~isscalar(linsys.C) && representsa_(params.V,'origin',eps) && ~any(params.vTrans)
    % y = Cx ... errors are scaled
    
    if isU
        Rout = zonotope(...
            linsys.C * (set.enc.c + set.boxFstartset_center + set.GuTrans_center), ...
            [linsys.C * set.enc.G, linsys.C * set.G_PUtotal_zero, ...
            linsys.C * diag(set.boxFstartset_Gbox + set.GuTrans_Gbox + set.G_PUtotal_infty)]);
    else
        Rout = zonotope(...
            linsys.C * (set.enc.c + set.boxFstartset_center + set.GuTrans_center), ...
            [linsys.C * set.enc.G, linsys.C * diag(set.boxFstartset_Gbox + set.GuTrans_Gbox)]);
    end
    Rout_error = errR2Y * Rcont_error;
    
    % time-point solution
    if isU
        Rout_tp = linsys.C * zonotope(set.Hstartp.c, ...
            [set.Hstartp.G, set.G_PUtotal_zero, diag(set.G_PUtotal_infty)]);
    else
        Rout_tp = linsys.C * set.Hstartp;
    end
    Rout_tp_error = errR2Y * Rcont_tp_error;

else
    % general case: y = Cx + V
    
    % compute output set and contained over-approximation error
    if isU
        Rout = zonotope(...
            linsys.C * (set.enc.c + set.boxFstartset_center + set.GuTrans_center) + params.vTrans, ...
            [linsys.C * generators(set.enc), linsys.C * set.G_PUtotal_zero, ...
            linsys.C * diag(set.G_PUtotal_infty + set.boxFstartset_Gbox + set.GuTrans_Gbox), ...
            generators(params.V)]);
    else
        Rout = zonotope(...
            linsys.C * (set.enc.c + set.boxFstartset_center + set.GuTrans_center) + params.vTrans, ...
            [linsys.C * set.enc.G, ...
            linsys.C * diag(set.boxFstartset_Gbox + set.GuTrans_Gbox), ...
            generators(params.V)]);
    end
    Rout_error = errR2Y * Rcont_error;
    
    % time-point solution
    if isU
        Rout_tp = zonotope(linsys.C * set.Hstartp.c + params.vTrans, ...
            [linsys.C * set.Hstartp.G, ...
            linsys.C * set.G_PUtotal_zero, ...
            linsys.C * diag(set.G_PUtotal_infty), generators(params.V)]);
    else
        Rout_tp = linsys.C * set.Hstartp + params.V + params.vTrans;
    end
    Rout_tp_error = errR2Y * Rcont_tp_error;
end

end


% error computation
function [errs,set,converged] = aux_errors(linsys,k,k_iter,errs,fullcomp,...
    isU,isuconst,G_U,u,set,t,timeStep,timeStepIdx,savedata)
% computes all errors for the given time step size, that is,
%   accumulating:     e.idv_PUtkplus1
%   non-accumulating: e.idv_linComb, e.idv_PUtauk, e.idv_F, e.idv_G
%   and reduction:    e.step_red (set to 0 as reduction is performed later)
% furthermore, auxiliary variables such as F, G are computed, as well
% as the set PU which is closely related to its respective error

% compute exponential matrix
eAdtk = computeField(linsys.taylor,'eAdt',struct('timeStep',timeStep));
eAtk = computeField(linsys.taylor,'eAdt',struct('timeStep',t));
converged = true;

% compute accumulating error
% - consists of eps_PU_tkplus1

% compute PU and eps_PU_tkplus1 (eps_PU_tauk only if fullcomp, see below)
if ~isU
    set.eAtkPU = zeros(linsys.nrOfStates,1);
    errs.idv_PUtkplus1{k,1}(k_iter,1) = 0;
elseif timeStepIdx == 0
    [converged,set.G_PU_zero{k_iter,1},set.G_PU_infty{k_iter,1},set.PU_A_sum_error{k_iter,1}] = ...
        aux_PU(linsys,G_U,timeStep);
    % if sums in computation did not converge, exit now
    if ~converged
        return;
    end
    % over-approximation of the Hausdorff distance between the exact
    % particular solution PU and the computed one
    errs.idv_PUtkplus1{k,1}(k_iter,1) = ...
        vecnorm( sum(abs( (eAtk * set.PU_A_sum_error{k_iter,1}) * G_U ),2) ) ...
        + vecnorm( sum(abs( eAtk * set.G_PU_infty{k_iter,1} ),2) );
    % note: G_PU_zero, G_PU_infty, PU_A_sum_error are saved outside
else
    % re-use values from memory (as time step size already used)
    set.G_PU_infty{k_iter,1} = savedata.G_PU_infty{timeStepIdx};
    set.PU_A_sum_error{k_iter,1} = savedata.PU_A_sum_error{timeStepIdx};
    errs.idv_PUtkplus1{k,1}(k_iter,1) = ...
        vecnorm( sum(abs( (eAtk * set.PU_A_sum_error{k_iter}) * G_U ),2) ) ...
        + vecnorm( sum(abs( eAtk * set.G_PU_infty{k_iter} ),2) );
    set.G_PU_zero{k_iter,1} = timeStep * G_U;
end


% compute non-accumulating error
% use NaN values by default
errs.idv_F{k,1}(k_iter,1) = NaN;
errs.idv_G{k,1}(k_iter,1) = NaN;
errs.idv_linComb{k,1}(k_iter,1) = NaN;
errs.idv_PUtauk{k,1}(k_iter,1) = NaN;
errs.seq_nonacc{k,1}(k_iter,1) = NaN;

if fullcomp
% - consists of eps_linComb + 2*err(F H*) + 2*err(G u) + eps_PU_tauk

    % linear combination error: eps_linComb
    errs.idv_linComb{k,1}(k_iter,1) = ...
        compute_eps_linComb(errs,eAdtk,set.startset);
    
    % curvature error from correction matrix for the state
    try
        F = correctionMatrixState(linsys,timeStep,Inf);
    catch ME
        if strcmp(ME.identifier,'CORA:notConverged')
            converged = false; return
        else
            rethrow(ME);
        end
    end
    [errs.idv_F{k,1}(k_iter,1), set.boxFstartset_center{k_iter,1}, set.boxFstartset_Gbox{k_iter,1}] = ...
        compute_eps_F(errs,F,set.startset);    
    
    % curvature error from correction matrix for the input
    set.GuTrans_center{k_iter,1} = zeros(linsys.nrOfStates,1);
    set.GuTrans_Gbox{k_iter,1} = zeros(linsys.nrOfStates,1);
    if ~any(u)
        errs.idv_G{k,1}(k_iter,1) = 0;
    else
        try
            G = correctionMatrixInput(linsys,timeStep,Inf);
        catch ME
            if strcmp(ME.identifier,'CORA:notConverged')
                converged = false; return
            else
                rethrow(ME);
            end
        end
        if isuconst && timeStepIdx > 0 && savedata.fullcomp(timeStepIdx)
            set.GuTrans_center{k_iter,1} = savedata.GuTrans_center{timeStepIdx};
            set.GuTrans_Gbox{k_iter,1} = savedata.GuTrans_Gbox{timeStepIdx};
            errs.idv_G{k,1}(k_iter,1) = savedata.eG{timeStepIdx};
        else
            [errs.idv_G{k,1}(k_iter,1), set.GuTrans_center{k_iter,1}, set.GuTrans_Gbox{k_iter,1}] = ...
                compute_eps_G(errs,G,u);
        end
    end
    
    
    % err(e^At PU)
    % over-approximation of the Hausdorff distance between the computed
    % particular solution PU and 0 (minimum solution over time interval)
    if ~isU
        errs.idv_PUtauk{k,1}(k_iter,1) = 0;
    else
        errs.idv_PUtauk{k,1}(k_iter,1) = ...
            vecnorm( sum(abs(eAtk * set.G_PU_zero{k_iter,1}),2) + ...
            sum(abs(eAtk * set.G_PU_infty{k_iter,1}),2) );
    end
    
    % total non-accumulating error
    errs.seq_nonacc{k,1}(k_iter,1) = errs.idv_linComb{k}(k_iter) ...
        + errs.idv_PUtauk{k}(k_iter) + errs.idv_F{k}(k_iter) + errs.idv_G{k}(k_iter);
end

% total accumulating error (contribution of current step)
% note that e.idv_PUtauk over-approximates e.idv_PUtkplus1, but for the
% accumulation to the next step, e.idv_PUtkplus1 is used
errs.step_acc{k,1}(k_iter,1) = errs.idv_PUtkplus1{k}(k_iter);

end

% particular solutions
function [conv,G_PU_zero,G_PU_infty,A_sum_error] = aux_PU(linsys,G_U,timeStep)
% computation of the particular solution due to the uncertain input set U
% as well as the contained over-approximation error with respect to the
% exact particular solution; using a Taylor series where the truncation
% order is increased until the additional values are so small that the
% stored number (finite precision!) does not change anymore
% we use only the generator matrix to spare calls to zonotope functions ->
% therefore, this function is copied from particularSolution_timeVarying

% initialize particular solution (only generator matrices)
G_PU_zero = timeStep * G_U;
PU_diag = sum(abs(G_PU_zero),2);

% initialize errors
G_PU_infty = [];
A_sum_error = zeros(linsys.nrOfStates);

eta = 1;
options = struct('timeStep',timeStep,'ithpower',eta);

% loop until floating-point precision
stoploop = false; conv = true;
while true

    options.ithpower = eta;
    Apower_mm = getTaylor(linsys,'Apower',options);
    options.ithpower = eta+1;
    dtoverfac = getTaylor(linsys,'dtoverfac',options);
    % additional term (only matrix)
    addTerm = Apower_mm * dtoverfac;
    addG_PU = addTerm * G_U;
    addG_PU_diag = sum(abs(addG_PU),2);
    
    % safety check (if time step size too large, then the sum converges
    % too late so we already have Inf values)
    % previous criterion: any(any(isinf(addPU_diag))) || any(any(isnan(addPU_diag)))
    % ... but very costly for large A!
    if eta == 75
        conv = false; return
    end
    
    % check if floating-point precision reached
    if all( abs(addG_PU_diag) <= eps * abs(PU_diag) )
        stoploop = true;
    end
    
    % add term to simplified value for convergence
    PU_diag = PU_diag + addG_PU_diag;
    
    % error terms
    A_sum_error = A_sum_error + addTerm;
    G_PU_infty = [G_PU_infty, addG_PU];
    
    % break loop
    if stoploop
        break;
    end
    
    % increment eta
    eta = eta + 1;

end

end

function [set,err_red_k] = aux_reduce(linsys,t,set,err_bound_red_k,isU)
% function which simultaneously propagates and reduces the particular
% solution due to the uncertainty U; we only work with generator matrices
% to spare function calls and object initializations

% overview of used variables
% G_eAtkPU_zero     generator matrix of e^Atk P^U_zero (for reduction)
% Gabs_eAtkPU_zero  absolute value of G_eAtkPU_zero
% G_eAtkPU_infty    generator matrix of e^Atk P^U_infty (-> box)
% G_PUtotal_zero    generator matrix of P^U_zero (for reduction)
% G_PUtotal_infty   generator matrix of P^U_infty (-> box)
% Gabs_PUtotal_zero absolute value of G_PUtotal_zero
% girard_PUtotal_zero   girard metric (norm_1 - norm_infty) of PUtotal_zero
% Ghat_PUtotal_zero generators of PUtotal_zero without maximum entry
%                   w.r.t infty-norm
% --- saved in set-struct for usage in successive steps
% ... G_PUtotal_zero nevers contains any girard-0-generators

if ~isU
    err_red_k = 0;
    return
end

% state dimension, propagation matrices
n = linsys.nrOfStates;
eAtk = readFieldForTimeStep(linsys.taylor,'eAdt',t);

% the _infty parts are always boxed (which induces no additional error),
% saved as a vector... correct handling in computation of output set
G_eAtkPU_infty = sum(abs( eAtk * set.G_PU_infty ),2);
set.G_PUtotal_infty = set.G_PUtotal_infty + G_eAtkPU_infty;

% generator matrix of additional solution e^Atk PU_zero
G_eAtkPU_zero = eAtk * set.G_PU_zero;
Gabs_eAtkPU_zero = abs(G_eAtkPU_zero);

% for fast falsification
set.GfastInner_add = linsys.C * G_eAtkPU_zero;

% reduction algorithm copied and adapted from zonotope/reduceAdaptive
% select generators using 'girard' metric (1-norm minus infty-norm)
girard_eAtkPU_zero = sum(Gabs_eAtkPU_zero,1) - max(Gabs_eAtkPU_zero,[],1);

% box axis-aligned generators
axisAligned = girard_eAtkPU_zero == 0;
if any(axisAligned)
    % add axis-aligned generators to infty part
    set.G_PUtotal_infty = set.G_PUtotal_infty ...
        + sum(Gabs_eAtkPU_zero(:,axisAligned),2);
    Gabs_eAtkPU_zero = Gabs_eAtkPU_zero(:,~axisAligned);
    G_eAtkPU_zero = G_eAtkPU_zero(:,~axisAligned);
end
nrG_eAtkPU_zero = size(G_eAtkPU_zero,2);

% compute auxiliary generators for new part G_eAtkPU_zero (without
% axis-aligned generators)
[maxval,maxidx] = max(Gabs_eAtkPU_zero,[],1);
% use linear indexing for speed-up
mu_eAtkPU_zero = zeros(n,nrG_eAtkPU_zero);
cols = n*(0:nrG_eAtkPU_zero-1);
mu_eAtkPU_zero(cols+maxidx) = maxval;
Ghat_eAtkPU_zero = Gabs_eAtkPU_zero - mu_eAtkPU_zero;

% reconstruct matrices for PUtotal_zero = PU_zero(tk) + e^Atk PU_zero(Delta tk)
set.G_PUtotal_zero(:,set.nrG_PUtotal_zero+1:set.nrG_PUtotal_zero+nrG_eAtkPU_zero) = G_eAtkPU_zero;
set.Ghat_PUtotal_zero(:,set.nrG_PUtotal_zero+1:set.nrG_PUtotal_zero+nrG_eAtkPU_zero) = Ghat_eAtkPU_zero;

set.girard_PUtotal_zero(:,set.nrG_PUtotal_zero+1:set.nrG_PUtotal_zero+nrG_eAtkPU_zero) = girard_eAtkPU_zero(~axisAligned);
[~,minidx] = min(set.girard_PUtotal_zero);

% compute over-approximation of dH (omega_max): only first value
% check whether this value is already larger than reduction error bound
omega_max_min = 2 * vecnorm(set.Ghat_PUtotal_zero(:,minidx),2);

if isempty(omega_max_min) || omega_max_min > err_bound_red_k
    % no generators can be reduced
    err_red_k = 0;
    set.nrG_PUtotal_zero = set.nrG_PUtotal_zero + nrG_eAtkPU_zero;
else
    set.nrG_PUtotal_zero = size(set.girard_PUtotal_zero,2);
    [~,idx] = mink(set.girard_PUtotal_zero,set.nrG_PUtotal_zero);
    
    % compute over-approximation of dH (omega_max)
    sum_temp = zeros(n,1);
    err_red_k = 0;
    redIdx = set.nrG_PUtotal_zero;
    for ijk=1:set.nrG_PUtotal_zero
        sum_temp = sum_temp + set.Ghat_PUtotal_zero(:,idx(ijk));
        omega_max_ijk = 2 * vecnorm(sum_temp,2);
        if omega_max_ijk > err_bound_red_k
            redIdx = ijk-1; break;
        else
            err_red_k = omega_max_ijk;
        end
    end
    
    % indices which are reduced
    idxRed = idx(1:redIdx);

    % add box made of generators selected for reduction to infty part
    set.G_PUtotal_infty = set.G_PUtotal_infty + sum(abs(set.G_PUtotal_zero(:,idxRed)),2);

    % remove indices which were selected for reduction
    set.girard_PUtotal_zero(:,idxRed) = [];
    set.Ghat_PUtotal_zero(:,idxRed) = [];
    set.G_PUtotal_zero(:,idxRed) = [];
    set.nrG_PUtotal_zero = set.nrG_PUtotal_zero - redIdx;
end

end

% adaptation of time step size
function maxTimeStep = aux_maxTimeStep(t,uTransVec,tu,tFinal,maxTimeStepSpec)
% maximum admissible value for time step size in the current step, which is
% determined by input vector u, which is required to be constant over
% the spanned time interval; additionally, the changing time of the
% specification satisfaction cannot be crossed in the time interval

% default: non time-varying input center
maxTimeStep = tFinal - t;

% time-varying input center
if size(uTransVec,2) > 1
    % add 1e-12 for numerical stability
    idxNextSwitch = find(tu > t + 1e-12,1,'first');
    % last step: idxNextSwitch = [], then remaining time as computed above
    if ~isempty(idxNextSwitch)
        maxTimeStep = tu(idxNextSwitch) - t;
    end
end

% change in specSAT boolean cannot be crossed
maxTimeStep = min([maxTimeStep, maxTimeStepSpec]);

end

function idx = aux_checkTimeStep(timeStep,savedata)
% checks whether auxiliary sets have already been computed for this time
% step size (allows to save a lot of re-computations)

idx = find(abs(savedata.timeStep - timeStep) < 1e-9,1);
if isempty(idx)
    idx = 0; % set 0 since this can be assigned to a list and is not an index
end

end

function [timeStep,timeStepIdx] = aux_adjustTimeStep(k,timeStep,isU,maxTimeStep,savedata)
% adjusts the predicted time step size in order to avoid recomputations;
% picks nearest smaller time step size from savedata, unless difference is
% larger than some fixed threshold ratio

if k == 1
    % different in init step... spend some time to find good time step size
    factor = 0.95;
elseif isU
    % limit decrease is quadratically over time step size
    % at least 50% of the error -> 75% of predicted time step size
    factor = 0.75;
else
    % limit decrease is linear -> achieve 80% of the error
    factor = 0.80;
end

% eligible time step sizes
idx_adm = savedata.timeStep < timeStep & savedata.timeStep > factor * timeStep;
% return if empty
if ~any(idx_adm)
    timeStepIdx = 0;
    % if proposed time step size between not largest or smallest until now,
    % set it to midway between two already existing ones
    if ~( all(timeStep >= savedata.timeStep) || all(timeStep <= savedata.timeStep) )
        timeStep_ub = min(savedata.timeStep(savedata.timeStep >= timeStep));
        timeStep_lb = max(savedata.timeStep(savedata.timeStep <= timeStep));
        timeStep = timeStep_lb + 0.5*(timeStep_ub - timeStep_lb);
    end
else
    % read time step and corresponding index from savedata
    [timeStep,temp] = max(savedata.timeStep(idx_adm));
    timeStepIdx = find(idx_adm,temp,'first');
    timeStepIdx = timeStepIdx(temp);
end

% ensure that maximum time step size not exceeded
if timeStep > maxTimeStep
    timeStep = maxTimeStep;
    timeStepIdx = max([0, find(abs(timeStep - savedata.timeStep) < 1e-14,1,'first')]);
end

end


% saving computations
function [timeStepIdx,savedata] = aux_savedata(savedata,...
    k,k_iter,fullcomp,isU,isuconst,set,errs,timeStep)
% save data into struct to avoid recomputations in subsequent steps
% this functions is entered in two cases:
% 1. the time step size has not been used before
% 2. the time step size has been used before, but fullcomp was false and
%    now fullcomp is true (so we have additional values to store)

timeStepIdx = find(abs(savedata.timeStep - timeStep) < 1e-14,1,'first');
if isempty(timeStepIdx)
    timeStepIdx = length(savedata.timeStep) + 1;
    savedata.Pu{timeStepIdx,1} = [];
end

% time step size
savedata.timeStep(timeStepIdx,1) = timeStep;
% affine dynamics
savedata.fullcomp(timeStepIdx,1) = fullcomp;
if isuconst
    if fullcomp
        savedata.GuTrans_center{timeStepIdx,1} = set.GuTrans_center{k_iter};
        savedata.GuTrans_Gbox{timeStepIdx,1} = set.GuTrans_Gbox{k_iter};
        savedata.eG{timeStepIdx,1} = errs.idv_G{k}(k_iter);
    else
        savedata.GuTrans_center{timeStepIdx,1} = [];
        savedata.GuTrans_Gbox{timeStepIdx,1} = [];
        savedata.eG{timeStepIdx,1} = [];
    end
end

% particular solution
if isU
    savedata.G_PU_zero{timeStepIdx,1} = set.G_PU_zero{k_iter};
    savedata.G_PU_infty{timeStepIdx,1} = set.G_PU_infty{k_iter};
    savedata.PU_A_sum_error{timeStepIdx,1} = set.PU_A_sum_error{k_iter};
else
    savedata.G_PU_zero{timeStepIdx,1} = [];
    savedata.G_PU_infty{timeStepIdx,1} = [];
    savedata.PU_A_sum_error{timeStepIdx,1} = [];
end

end

% ------------------------------ END OF CODE ------------------------------
