function [timeInt,timePoint,res,savedata] = reach_adaptive(obj,options)
% reach_adaptive - computes an outer-approximation of the reachable set
%    for linear time-invariant systems given a maximum error in terms of
%    the Hausdorff distance to the exact reachable set; all internal
%    reachability settings are set automatically
%
% Syntax:
%    [Rout,Rout_tp,res,savedata] = reach_adaptive(obj,options)
%
% Inputs:
%    obj - linearSys object
%    options - model parameters and options for the computation of reachable sets
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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       08-July-2021
% Last update:   22-March-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initializations ---------------------------------------------------------

% show debug help (only dev)
isdebug = false;

% call from verify?
if ~isfield(options,'verify')
    options.verify = false;
end

% init step and time (shift so that start at 0), and spec satisfaction
k = 0; t = 0; res = true;
options.tFinal = options.tFinal - options.tStart;

% init struct which keeps sets/expmat/timeStep in order to avoid
% recomputations (possibly also useful as a return argument)
savedata = aux_initSavedata(options);

% time intervals where specifications are not yet verified
% TODO: integrate in validateOptions later...
[options.specUnsat,options.tFinal] = aux_initSpecUnsat(options,savedata);

% initialize time step vector and previous time step (dummy value)
tVec = -1;

% initialize all variables for sets
[set,Rcont_error,Rcont_tp_error,timeInt,timePoint] = aux_initSets(obj,options);

% initialize all variables for exponential (propagation) matrix
expmat = aux_initExpmat(obj);

% rewrite equations in canonical form
options = aux_canonicalForm(obj,options);

% init all auxiliary structs and flags
[e,ebar,isU,G_U,isuconst,coeff,bisect,quickCheck,options] = ...
    aux_initStructsFlags(obj,options);

% pre-compute allocation of reduction error
[ebar,savedata] = aux_ebarred(obj,isU,G_U,e,ebar,options,savedata);

% struct to store near-all values (helpful for debugging)
debugdata = [];

% compute first output solution and corresponding error
timePoint.set{1} = obj.C * options.R0 + options.V + options.vTransVec(:,1);
timePoint.error(1) = 0;
% -------------------------------------------------------------------------


% main loop until the time is up
while options.tFinal - t > 1e-9
    
    % update step counter
    k = k + 1;
    
    % log information
    verboseLog(k,t,options);
    
    % update constant input vector (if necessary)
    [options.uTrans, options.vTrans] = aux_uTrans_vTrans(t,...
        options.uTransVec,options.vTransVec,options.tu);
    
    % flag whether full computation required (only where specifications are 
    % not yet satisfied), additionally return time until next change
    [fullcomp(k,1), maxTimeStepSpec] = aux_fullcomp(t,options.specUnsat);
    
    % compute maximum admissible time step size for current step (either
    % until the time horizon or the delayed start, or shorter if there is a
    % time-varying input vector)
    maxTimeStep = aux_maxTimeStep(t,options.uTransVec,options.tu,...
        options.tFinal,maxTimeStepSpec);
    
    % take care of initializations in first step
    if k == 1
        % initialize start set
        set.startset = options.R0;
        
        % initialize expmat
        expmat.tk = eye(obj.dim);

        % guess initial time step size as the time horizon
        timeStep = min([maxTimeStep, options.tFinal]);
        
        % time step adaptation using approximation functions or bisection:
        % first step always bisection (more robust)
        approxfuncs = false;
        
    else
        if fullcomp(k)
        % define startset for current step
        set.startset = set.Hstartp;
        end
        
        % reuse same time step size as in the previous step (unless this
        % exceeds the value for the maximum admissible timeStep)
        timeStep = min([maxTimeStep, timeStep]);
        
        % time step adaptation using approximation functions or bisection:
        approxfuncs = true;
    end
    
    % initialize current step struct
    cnti = 0;
    cnt = [];
    
    % loop until both errors are satisfied
    while true
        
        % increase counter
        cnti = cnti + 1;
        
        % compute error bounds
        ebar = aux_ebar(t,k,e,ebar,timeStep,options.tFinal);
        
        % check if time step size was used before: note that this call is
        % not redundant if maxTimeStep used or in first step
        cnt.timeStepIdx(cnti,1) = aux_checkTimeStep(timeStep,savedata);
        timeStepIdxs(k,1) = cnti;
        
        % check if time step size is the same as in the last step
        if k > 1 && cnti == 1 && isU
            expmat.tk = expmat.tkplus1;
        end
        
        % compute error bounds and errors for given timeStep
        [e,set,expmat] = aux_errors(fullcomp(k),obj,isU,isuconst,G_U,...
            options.uTrans,set,e,expmat,timeStep,cnt.timeStepIdx(cnti,1),savedata);
        % if time step size is too big, then Inf in aux_F_Ftilde
        % -> repeat with smaller time step size
        if ~expmat.conv
            if cnti == 1
                timeStep = timeStep * 0.1; cnti = 0; continue
            else
                % take previous time step (has to be ok)
                [~,set,e,ebar,expmat,timeStep] = aux_ReadWriteLoop(k,cnt,...
                    fullcomp(k),isU,set,e,ebar,expmat,timeStep,cnti-1,'read');
                timeStepIdxs(k,1) = find(cnt.timeStepIdx(1:cnti) == cnt.timeStepIdx(cnti-1),1,'first');
                break
            end
        else
            % save results
            cnt = aux_ReadWriteLoop(k,cnt,fullcomp(k),isU,...
                set,e,ebar,expmat,timeStep,cnti,'write');
            % if new time step size used, save auxiliary variables
            if cnt.timeStepIdx(cnti) == 0 || fullcomp(k) ~= savedata.fullcomp(cnt.timeStepIdx(cnti))
                [cnt.timeStepIdx(cnti,1),savedata] = aux_savedata(...
                    fullcomp(k),isU,isuconst,set,e,expmat,timeStep,savedata);
            end
        end

        % init/update coefficients of approximation functions
        if approxfuncs
            coeff = aux_coeff(cnti,fullcomp(k),coeff,e,timeStep,cnt);
        end
        
        % debug: print fulfillment of error bounds
        aux_print_e2ebar(isdebug,k,cnti,isU,fullcomp(k),e,ebar,timeStep);
        
        % debug: save data (required for bounds check and plot functions)
        debugdata = aux_debug_debugdata(k,cnti,fullcomp(k),debugdata,e,ebar,coeff,timeStep);
        
        % check non-accumulating and accumulating error against bound
        cnt.accok(cnti,1) = e.acc < ebar.acc(k);
        cnt.nonaccok(cnti,1) = e.nonacc < ebar.rem(k) || ~fullcomp(k);
        
        % update bisect values (lower/upper bounds also used in adaptation
        % approximation using functions)
        bisect = aux_bisect(k,cnti,bisect,isU,e,ebar,timeStep,cnt);
        
        % check whether suitable time step size found
        if cnt.accok(cnti,1) && (~fullcomp(k) || cnt.nonaccok(cnti,1))
            
            % proceed if timeStep is already maximum admissible time step
            if abs(timeStep - maxTimeStep) < 1e-12
                break;
            
            % ensure that either the joint bound for eacc and enonacc
            % or the bound for eacc is fulfilled relatively closely
            elseif e.acc > 0.90 * ebar.acc(k) || (fullcomp(k) && ...
                    e.acc + e.nonacc - e.PU_tkplus1 > 0.90 * ebar.rem(k) )
                break;
                
            end
        end

        % only for aux_coeff
        cnt.ePU_tkplus1(cnti,1) = e.PU_tkplus1;
        cnt.eF(cnti,1) = e.F;
        cnt.eFtilde(cnti,1) = e.Ftilde;
        
        % predict timeStep necessary to satisfy both (non-acc/acc) bounds
        if approxfuncs
            timeStep = aux_predTimeStep_approxfuncs(t,k,fullcomp(k),coeff,...
                e,ebar,timeStep,maxTimeStep,options.tFinal);
            % approx. function method proposes a value for the time step
            % size which we cannot accept if:
            % 1) value smaller than any lower bound (where errors fulfilled)
            % 2) value larger than any upper bound (where errors not fulfilled)
            % 3) value smaller than any upper bound (where errors fulfilled)
            approxfuncs = ~( ...
                any(timeStep <= [bisect.lb.timeStep_acc, bisect.lb.timeStep_nonacc]) || ...
                any(timeStep >= [bisect.ub.timeStep_acc, bisect.ub.timeStep_nonacc] & ...
                    ~[bisect.ub.accok, bisect.ub.nonaccok]) || ...
                any(timeStep <= [bisect.ub.timeStep_acc, bisect.ub.timeStep_nonacc] & ...
                    [bisect.ub.accok, bisect.ub.nonaccok]) );
        end
        
        if ~approxfuncs
            timeStep = aux_predTimeStep_bisect(t,k,fullcomp(k),bisect,isU,...
                e,ebar,timeStep,maxTimeStep,options.tFinal);
        end
        
        % adjust chosen time step size to avoid recomputations
        [timeStep,cnt.timeStepIdx(cnti+1,1)] = ...
            aux_adjustTimeStep(k,timeStep,isU,maxTimeStep,savedata);
        
        % exit if same time step size proposed as computed in a prior cnti
        timeStepIdx_ = find(cnt.timeStepIdx(1:cnti) == cnt.timeStepIdx(cnti+1),1,'first');
        if ~isempty(timeStepIdx_) && cnt.accok(timeStepIdx_) ...
                && (~fullcomp(k) || cnt.nonaccok(timeStepIdx_))
            timeStepIdxs(k,1) = timeStepIdx_;
            [~,set,e,ebar,expmat,timeStep] = aux_ReadWriteLoop(k,cnt,...
                fullcomp(k),isU,set,e,ebar,expmat,timeStep,timeStepIdx_,'read');
            break;
        end
        
    end

    % compute propagation matrix
    if isU
        expmat.tkplus1 = expmat.Deltatk * expmat.tk;
    elseif ~fullcomp(k) && abs(timeStep - maxTimeStepSpec) < 1e-9
        % switch from partial computation to full computation in next step
        expmat.tkplus1 = expm(obj.A*(t+timeStep));
    end
    
    % accumulate particular solution
    e.red = 0;
    if isU
        [set,e] = aux_reduce(obj,k,set,e,ebar,expmat); % omega_max
%         [set,e] = aux_reduce_rad(obj,k,set,e,ebar,expmat); % omega_rad
    end
    
    % accumulate accumulating and reduction error
    if k == 1
        e.acctotal(k,1) = e.acc;
        e.redtotal(k,1) = e.red;
    else
        e.acctotal(k,1) = e.acctotal(k-1) + e.acc;
        e.redtotal(k,1) = e.redtotal(k-1) + e.red;
    end
    
    
    if any(any(options.uTransVec))
    % compute particular solution due to input vector
    % due to current method of computation, this induces no error as the
    % solution is just a point!
    if ~isuconst || isempty(savedata.Pu{cnt.timeStepIdx(timeStepIdxs(k))})
        [set.Pu,expmat] = aux_Pu(obj,options.uTrans,expmat,timeStep);
        if isuconst
            savedata.Pu{cnt.timeStepIdx(timeStepIdxs(k)),1} = set.Pu;
        end
    else
        % spare re-computation if time step size used before
        set.Pu = savedata.Pu{cnt.timeStepIdx(timeStepIdxs(k))};
    end
    
    % propagate constant input solution
    if k == 1
        set.Putotal = set.Pu;
    else
        set.Putotal = expmat.Deltatk * set.Putotal + set.Pu;
    end
      
    end
    
    % propagate reachable sets
    if ~fullcomp(k)
        % no computation of the full affine solution or reachable sets/output
        % sets (incl. the contained errors)
        
        % save NaN value for time-interval errors (to ease debugging)
        Rcont_error(k,1) = NaN;
        
        % compute time-point affine solution at the end of the step, which is
        % required to start the full computation in the next step
        if abs(timeStep - maxTimeStepSpec) < 1e-9
            set.Hstartp = expmat.tkplus1*options.R0 + set.Putotal;
        end
        
        % set Rout entry to empty and error to NaN for consistency with tVec
        timeInt.set{k,1} = [];
        timeInt.error(k,1) = NaN;
        
        % propagate auxiliary term for quick check of inner-approximations
        if quickCheck
            [res,set,savedata] = aux_quickCheck(obj,set,[t,t+timeStep],[],...
                e.redtotal(k),Rout_error(k),options,savedata);
        end
    
    else
        % full propagation of sets and errors
            
        % compute auxiliary solution of affine system
        set.Hstartp = expmat.Deltatk * set.startset + set.Pu;
        set.enc = enclose(set.startset, set.Hstartp);
        
        % gather over-approximation error (time-point/interval)
        if k == 1
            Rcont_error(k,1) = e.nonacc + e.redtotal(k);
        else
            Rcont_error(k,1) = e.nonacc + e.acctotal(k-1) + e.redtotal(k);
        end
        Rcont_tp_error(k,1) = e.acctotal(k) + e.redtotal(k);

        % compute output set and contained over-approximation error
        [timeInt.set{k,1},timeInt.error(k,1),timePoint.set{k+1,1},timePoint.error(k+1,1),set] = ...
            aux_outputset(obj,isU,set,Rcont_error(k),options,Rcont_tp_error(k));
        if options.verify && quickCheck
            % check violation of inner-approximation for quick exit
            [res,set,savedata] = aux_quickCheck(obj,set,[t,t+timeStep],...
                center(timeInt.set{k}),e.redtotal(k),timeInt.error(k),options,savedata);
        end
    
    end
    
    % update debugdata
    debugdata = aux_debug_debugdata_sets(k,fullcomp,debugdata,set,expmat);
    
    % save chosen time step to time step vector
    tVec(k,1) = timeStep;
    
    % update starting time of next step (needs to be done before checking
    % the loop condition with respect to the time horizon)
    t = t + timeStep;
    
    % exit if simple-to-check specifications are violated
    if ~res
        break
    end
    
end

% log information
verboseLog(k,t,options);

% debug: plot errors and bounds
% aux_plot_errorcurve(isdebug,fullcomp,e,ebar,tVec,debugdata);

% debug: plot fulfillment of error bounds
% aux_plot_e2ebar(isdebug,fullcomp,debugdata,timeStepIdxs);

% debug: plot error curve
% aux_plot_Rerror(isdebug,fullcomp,tVec,Rcont_error,Rcont_tp_error,e,options.tu);

% check if all errors below respective bounds
% aux_check_errors(fullcomp,e,ebar,Rcont_error,Rcont_tp_error,tVec,debugdata,timeStepIdxs);

% write time to output structs
timePoint.time = num2cell([0;cumsum(tVec)]+options.tStart);
timeInt.time = cell(length(timePoint.time)-1,1);
for k=1:length(timePoint.time)-1
    timeInt.time{k} = interval(timePoint.time{k},timePoint.time{k+1});
end

% specification fulfilled (currently, no specifications included)
res = true;

end


% Auxiliary functions -----------------------------------------------------

% main function to check correct functionality
function aux_check_errors(fullstep,e,ebar,Rout_error,Rout_tp_error,tVec,debugdata,timeStepIdxs)
% check function which returns true/false depending on fulfillment of error
% bounds over time; as a safety measure, this function is always executed

fullstep_tp = fullstep;
fullstep_tp(find(~fullstep_tp,1,'last')) = true;

% assume checks to be ok
checkok = true;
internalcheck = true;

steps = length(debugdata.e.nonacc);
% bound of non-accumulating error
ebarnonacc = ebar.rem;
enonacc = zeros(steps,1); eacc = zeros(steps,1); ered = zeros(steps,1);
for i=1:steps
    enonacc(i,1) = debugdata.e.nonacc{i}(timeStepIdxs(i));
    eacc(i,1) = debugdata.e.acc{i}(timeStepIdxs(i));
    ered(i,1) = debugdata.e.red{i}(timeStepIdxs(i));
end
% quick fix for time-point error (empty in verify-call)
if isempty(Rout_tp_error)
    Rout_tp_error = zeros(steps,1);
end

% cumulative sum to yield bound for accumulating error
% note: subtract error committed in the step!
ebaracctotal = e.acctotal + ebar.acc - eacc;
% ebarredtotal = e.redtotal + ebarred - ered;

% compute total committed error
totalerror = enonacc + e.acctotal - eacc + e.redtotal;
totalerror_tp = e.acctotal + e.redtotal;

% --- checks ---
% all errors and bounds need to be larger or equal to zero
% (just to find obvious bugs)
if any([Rout_error(fullstep); Rout_tp_error(fullstep_tp); enonacc(fullstep); ...
        eacc; ered; ebarnonacc; ebar.acc; ebar.red; ebaracctotal] < 0)
    checkok = false;
end

% full errors (time-interval and time-point)
% (1) need to be below maximum error
if any(Rout_error(fullstep) > e.emax) || any(Rout_tp_error(fullstep_tp) > e.emax)
    checkok = false;
end
% (2) re-computation of full errors has to match ongoing computation
%     note: use a tolerance here...
if any(abs(totalerror(fullstep) - Rout_error(fullstep)) > 1e-9) || ...
        (any(Rout_tp_error) && any(abs(totalerror_tp(fullstep_tp) - Rout_tp_error(fullstep_tp)) > 1e-9))
    internalcheck = false;
end

% non-accumulating errors (linComb, F, Ftilde)
% (1) need to be below maximum error
% (2) need to be below remaining error (after subtraction of acc error
%     until now and reduction error bound)
if any(enonacc > e.emax)
    checkok = false;
end
if any(enonacc > ebarnonacc)
    internalcheck = false;
end

% accumulating error (PU)
% (1) single step error needs to be below accumulated error
% (2) single step error needs to be below corresponding bound
% (3) accumulated error needs to be below maximum error
% (4) accumulated error needs to be below final admissible value
if any(e.acctotal > e.emax) || any(e.acctotal > e.emax - ebar.red_e(end))
    checkok = false;
end
if any(eacc > e.acctotal) || any(eacc > ebar.acc) 
    internalcheck = false;
end
% (5) accumulated error needs to be below linearly increasing bound
%     (note: only in the current version!)
% (6) accumulated error bound needs to be below linearly increasing bound
% compute bound first...
finalValue = e.emax - ebar.red_e(end);
cumsumtVec_tp = cumsum(tVec);
linBound = finalValue * cumsumtVec_tp / cumsumtVec_tp(end);
if any(e.acctotal > linBound) || any(linBound - ebaracctotal < -10*eps)
    internalcheck = false;
end

% reduction error (PU)
% (1) single step error needs to be below accumulated error
% (2) single step error needs to be below corresponding bound
% (3) accumulated error needs to be below maximum error
% (4) accumulated error needs to be below final admissible value
if any(e.redtotal > e.emax)
    checkok = false;
end
if any(ered > e.redtotal) || any(ered > ebar.red) || ...
         any(e.redtotal > ebar.red_e(end))
    internalcheck = false;
end


% print result of checks to console
disp("--------------------------------");
% is maximum error bound respected of the whole time horizon?
if checkok
disp("   Error bounds respected");
else
disp("   Error bounds NOT respected");
end
% do the internal errors satisfy their bounds?
if internalcheck
disp("   (Internal checks ok)");
else
disp("   (Internal checks failed)");
end
disp("--------------------------------");

end

% initializations (all before main loop)
function [set,Rcont_error,Rcont_tp_error,timeInt,timePoint] = aux_initSets(obj,options)
% init over-approximation error w.r.t the corresponding exact reachable
% sets (this is not a return value, but might still be convenient to have)
Rcont_error = [];
Rcont_tp_error = [];

% initialize time-point and time-interval output sets and
% over-approximation error w.r.t the corresponding exact output sets
% remark: index of Rout_tp: 1 ... t=t_0, 2 ... t=t_1, etc.
%         index of Rout:    1 ... t=[t_0,t_1], 2 ... [t_1,t_2]
% ... hence, Rout_tp has one more entry compared to Rout
timeInt.set = cell(0);
timeInt.time = cell(0);
timeInt.error = [];
timePoint.set = cell(0);
timePoint.time = cell(0);
timePoint.error = [];

% initialize auxiliary solutions
set.startset = [];
set.Hstarti = [];
set.Hstartp = [];
set.Pu = zeros(obj.dim,1);
set.Putotal = zeros(obj.dim,1);
set.PU = [];
set.PUtotal = [];
set.Rcont_tp_kminus1 = options.R0;
% for reduction and propagation of particular solution
set.girard_PUtotal_zero = [];
set.G_PUtotal_infty = zeros(obj.dim,1);
set.G_PUtotal_zero = double.empty(obj.dim,0);
set.Gabs_PUtotal_zero = double.empty(obj.dim,0); % only omega_rad
set.Ghat_PUtotal_zero = double.empty(obj.dim,0); % only omega_max
set.nrG_PUtotal_zero = 0;
% set.nrG_PUtotal_zero_incr = floor(1024^2/8/obj.dim); % ~1MB storage space
% set.nrG_PUtotal_zero_full = set.nrG_PUtotal_zero_incr;
% set.G_PUtotal_zero = zeros(obj.dim,set.nrG_PUtotal_zero_full);
% set.Gabs_PUtotal_zero = zeros(obj.dim,set.nrG_PUtotal_zero_full); % only omega_rad
% set.Ghat_PUtotal_zero = zeros(obj.dim,set.nrG_PUtotal_zero_full); % only omega_max
% only fast inner-approximation check
if isfield(options,'savedata')
    if isfield(options.savedata,'Gunsat')
    for i=1:length(options.savedata.Gunsat)
        if options.G{i}.fastInner
            set.G_GfastInner{i} = 0;
        end
    end
    end
    if isfield(options.savedata,'Funsat')
    for i=1:length(options.savedata.Funsat)
        if options.F{i}.fastInner
            set.F_GfastInner{i} = 0;
        end
    end
    end
end
set.GfastInner_add = zeros(obj.nrOfOutputs,1);

end

function options = aux_canonicalForm(obj,options)
% put inhomogeneity to canonical forms:
%    Ax + Bu + c + w  ->  Ax + u, where u \in U + uTransVec
%    Cx + Du + k + v  ->  Cx + v, where v \in V + vTransVec

% read-out disturbance and sensor noise
centerW = center(options.W);
W = options.W + (-centerW);
centerV = center(options.V);
V = options.V + (-centerV);

% initialize input vector if sequence given
if isfield(options,'uTransVec')
    % time-varying input vector
    uVec = options.uTransVec;
else
    % no time-varying input vector, but uTrans given
    uVec = options.uTrans;
end

% output equation
if any(any(obj.D))
    options.V = obj.D * options.U + V;
else
    options.V = V;
end

% state equation
options.U = obj.B * options.U + W;
options.uTransVec = obj.B * uVec + obj.c + centerW;

% initialize input vector for output equation
% note: this CANNOT the same as the input vector for the state equation
options.vTransVec = obj.k + centerV;

% assign first values
options.uTrans = options.uTransVec(:,1);
options.vTrans = options.vTransVec(:,1);

% remove fields for safety
options = rmfield(options,'W');
% note: U and V now overwritten!

end

function expmat = aux_initExpmat(obj)

% init convergence for exponential matrix auxiliary terms
expmat.conv = true;
% initialize power of A and |A|
expmat.Apower{1} = obj.A;
expmat.Apower_abs{1} = abs(obj.A);
% initialize positive and negative indices of A^eta
expmat.Apos = cell(0);
expmat.Aneg = cell(0);
% initialize A-related values: propagation matrices, powers, inverse
expmat.Deltatk = [];     % propagation from start to end of current step
expmat.tk = [];          % propagation until start of current step
expmat.tkplus1 = [];     % propagation until end of current step

% precompute inverse of A matrix
expmat.isAinv = rank(full(obj.A)) == obj.dim;
expmat.Ainv = [];
if expmat.isAinv
    expmat.Ainv = inv(obj.A);
end

end

function [e,ebar,isU,G_U,isuconst,coeff,bisect,quickCheck,options] = aux_initStructsFlags(obj,options)

% saving of operations for affine systems (u = const. over entire time
% horizon) vs. system with varying u or even uncertainty U
% -> the resulting if-else branching looks quite ugly, but still yields
% large speed-ups for high-dimensional systems
G_U = generators(options.U);
isU = ~isempty(G_U);
isuconst = ~(isfield(options,'uTransVec') && size(options.uTransVec,2) > 1);
% isaffineuconst = size(options.uTransVec,2) == 1 && ~isU;
% sparsity for speed up (acc. to numeric tests only for very sparse
% matrices actually effective)
if nnz(obj.A) / numel(obj.A) < 0.1
    obj.A = sparse(obj.A);
end

% compute factor for scaling of Hausdorff distance from state space to
% output space and vice versa due to linear transformation with C
options.errR2Y = norm(sqrtm(full(obj.C * obj.C')),2); % same as norm(obj.C)
% scale options.error so that it corresponds to R (this is only relevant in
% case an output equation is provided as options.errR2Y = 1 otherwise)
options.error = options.error / options.errR2Y;

% reassign maximum error
e.emax = options.error;

% initialize errors (\varepsilon)
e.nonacc = [];
e.acc = [];
e.red = [];
e.acctotal = [];
e.redtotal = [];
% e.norm2 = []; e.normF = []; e.normsqrt = [];

% initialize error bounds (\bar{\varepsilon}})
ebar.nonacc = [];
ebar.acc = [];
ebar.red = [];

% initialize coefficients of approximation function
coeff = [];

% initialize bisect struct
bisect = [];

% quick check of specs
quickCheck = false;
if isfield(options,'G') && ~isempty(options.G)
   Gtmp = cell2mat(options.G); 
   if any([Gtmp.fastInner])
        quickCheck = true;
   end
end
if isfield(options,'F') && ~isempty(options.F)
   Ftmp = cell2mat(options.F);
   if any([Ftmp.fastInner])
        quickCheck = true;
   end
end

end

function [ebar,savedata] = aux_ebarred(obj,isU,G_U,e,ebar,options,savedata)
% heuristics to determine a curve (params.tStart always shifted to 0!)
%   x-axis (time):  0 -> params.tFinal
%   y-axis (error): 0 -> (maximum) e.emax
% for the reduction error to yield the smallest zonotope order at the end
% of the time horizon; note that we only require to reduce the particular
% solution due to the uncertainty in the input (options.U); the result is a
% curve from which we can read the admissible reduction error and thus also
% the remaining error margin for the other errors for any point in time 
% (this function is only called once after the initial step)

if ~isU
    % no reduction, entire error margin is available for nonacc errors
    ebar.red_t = [0;options.tFinal];
    ebar.red_e = [0;0];
    savedata.reductionerrorpercentage = ebar.red_e(end)/e.emax;
    return
elseif isfield(options,'reductionerrorpercentage')
    % define error curve manually
    ebar.red_t = [0;options.tFinal];
    ebar.red_e = [0;e.emax * options.reductionerrorpercentage];
    savedata.reductionerrorpercentage = ebar.red_e(end)/e.emax;
    return
elseif isfield(savedata,'reductionerrorpercentage')
    % reduction error defined by previous run (only verify)
    ebar.red_t = [0;options.tFinal];
    ebar.red_e = [0;e.emax * options.savedata.reductionerrorpercentage];
    return
end
% hard-coded allocations:
% ebar.red_e = [0;0];
% half of entire error margin
% ebar.red_e = [0;e.emax/2];
% post-step: vector over time for reduction error allocation

% heuristics for near-optimal allocation
stepsforebar = 100;
timeStep = options.tFinal / stepsforebar;
expmat_aux = eye(obj.dim);
expmat_Deltatk = expm(obj.A * timeStep);

% 1. compute auxiliary sets V_k and errors e(V_k)
V = cell(stepsforebar,1);
errV = zeros(stepsforebar,1);

DeltatkU = timeStep * G_U;
for i=1:stepsforebar
    % compute auxiliary set
    V{i} = expmat_aux * DeltatkU;
    % propagate exponential matrix for next V
    expmat_aux = expmat_aux * expmat_Deltatk;
    % compute error (center is always zero)
    errV(i) = vecnorm( sum(abs(V{i}),2) );
end

% 2. compute weights and ordering
weights = errV ./ sum(errV);
% TODO: clarify how are ties decided?
[~,tau] = sort(weights);
% 3. sort V
errVsort_cumsum = cumsum(errV(tau));

% question: how much could we reduce if we allocate some portion of emax?
meshsize = 1000;
emax_percentage4ered = linspace(0,1,meshsize)';
emax_percentage4ered = emax_percentage4ered(1:end-1);
finalorder = zeros(meshsize-1,1);
% best final order is a priori the one where we do not reduce at all
bestfinalorder = stepsforebar+1;
min_idx = 1;

% loop over the mesh of emax to find a near-optimal result
for i=1:meshsize-1
    % find chi* = number of V_k that can be reduced for current red error
    idx = find(errVsort_cumsum < e.emax * emax_percentage4ered(i),1,'last');
    % fraction emax_percentage4ered is allocated to the reduction error,
    % results in factor N of total number of steps
    N = 1/(1-emax_percentage4ered(i));
    % compute resulting zonotope order
    if isempty(idx)
        finalorder(i) = N*(stepsforebar+1);
    else
        finalorder(i) = N*(stepsforebar+1 - idx);
    end
    % best reduction error allocation is the one which results in the
    % lowest zonotope order, also save the corresponding idx
    if finalorder(i) < bestfinalorder
        bestfinalorder = finalorder(i);
        min_idx = idx;
    end
end

% use sets V_k to generate a curve of admissible reduction errors over time
% ebar.red_t = (timeStep:timeStep:tFinal)';
% % append time horizon at the end of time vector
% if abs(ebar.red_t(end) - tFinal) > 10*eps
%     ebar.red_t = [ebar.red_t; tFinal];
% end
% ebar.red_e = zeros(length(ebar.red_t),1);
% for i=1:min_idx
%     ebar.red_e(tau(i)) = errV(tau(i));
% end
% ebar.red_e = cumsum(ebar.red_e);
% add zero at the start
% ebar.red_t = [0;ebar.red_t];
% ebar.red_e = [0;ebar.red_e];

% simpler method
ebar.red_t = [0;options.tFinal];
ebar.red_e = [0;e.emax * emax_percentage4ered(min_idx)];

% save value
savedata.reductionerrorpercentage = ebar.red_e(end)/e.emax;

% previous method... (U = options.U)
% % note: omit last step for now...
% % 1. compute auxiliary sets V_k and errors e(V_k)
% steps = 1000; % TODO: adapt this value
% timeStep = tFinal / steps;
% expmat_aux = eye(obj.dim);
% expmat_Deltatk = expm(obj.A * timeStep);
% 
% V = cell(steps,1);
% errV = zeros(steps,1);
% 
% DeltatkU = timeStep * U;
% for i=1:steps
%     % compute auxiliary set
%     V{i} = expmat_aux * DeltatkU;
%     % propagate exponential matrix for next V
%     expmat_aux = expmat_aux * expmat_Deltatk;
%     % compute error
%     errV(i) = aux_errOp(V{i});
% end
% 
% % 2. compute weights and ordering
% weights = errV ./ sum(errV);
% % TODO: clarify how are ties decided?
% [~,tau] = sort(weights);
% % 3. sort V
% errVsort_cumsum = cumsum(errV(tau));
% 
% % alternative (but practically not feasible)
% % Vsort = V(tau);
% % errVsort_cumsum = zeros(steps,1);
% % errVsort_cumsum(1) = aux_errOp(Vsort{1});
% % V_sum_temp = Vsort{1};
% % for i=2:steps
% %     V_sum_temp = V_sum_temp + Vsort{2};
% %     errVsort_cumsum(i) = aux_errOp(V_sum_temp);
% % end
% 
% % question: how much could we reduce if we allocate some portion of emax?
% meshsize = 1000;
% emax_cont = linspace(0,e.emax,meshsize)';
% finalorder = zeros(meshsize,1);
% % best final order is a priori the one where we do not reduce at all
% bestfinalorder = steps+1;
% 
% % loop over the mesh of emax to find a near-optimal result
% for i=1:meshsize-1
%     % find chi* = number of V_k that can be reduced for current red error
%     idx = find(errVsort_cumsum < emax_cont(i),1,'last');
%     % fraction 1/N is allocated to the reduction error
%     N = 1/(1-emax_cont(i));
%     % compute resulting zonotope order
%     if isempty(idx)
%         finalorder(i) = N*(steps+1);
%     else
%         finalorder(i) = N*(steps+1) - idx;
%     end
%     % best reduction error allocation is the one which results in the lowest
%     % zonotope order, also save the corresponding idx
%     if finalorder(i) < bestfinalorder
%         bestfinalorder = finalorder(i);
%         min_idx = idx;
%     end
% end
% 
% % use sets V_k to generate a curve of admissible reduction errors over time
% ebar.red_t = (timeStep:timeStep:tFinal)';
% % append time horizon at the end of time vector
% if abs(ebar.red_t(end) - tFinal) > 10*eps
%     ebar.red_t = [ebar.red_t; tFinal];
% end
% ebar.red_e = zeros(length(ebar.red_t),1);
% for i=1:min_idx
%     ebar.red_e(tau(i)) = errV(tau(i));
% end
% ebar.red_e = cumsum(ebar.red_e);
% % add zero at the start
% ebar.red_t = [0;ebar.red_t];
% ebar.red_e = [0;ebar.red_e];

end

function savedata = aux_initSavedata(options)
% initialize all fields in struct savedata
% sets/expmat correspond to timeStep of same index

% savedata already there from previous run
if isfield(options,'savedata')
    savedata = options.savedata;
else
    savedata = struct();
end

if ~isfield(savedata,'timeStep')
    % time step size
    savedata.timeStep = [];
end
if ~isfield(savedata,'expmatDeltatk')
    % propagation matrix
    savedata.expmatDeltatk = {};
end

% affine dynamics
if ~isfield(savedata,'fullcomp')
    savedata.fullcomp = [];
end
if ~isfield(savedata,'F')
    savedata.F = {};
end
if ~isfield(savedata,'Frad')
    savedata.Frad = {};
end
if ~isfield(savedata,'Fcenter')
    savedata.Fcenter = {};
end
if ~isfield(savedata,'Ftilde')
    savedata.Ftilde = {};
end
if ~isfield(savedata,'Ftilderad')
    savedata.Ftilderad = {};
end
if ~isfield(savedata,'Ftildecenter')
    savedata.Ftildecenter = {};
end
if ~isfield(savedata,'FtildeuTrans_center')
    savedata.FtildeuTrans_center = {};
end
if ~isfield(savedata,'FtildeuTrans_Gbox')
    savedata.FtildeuTrans_Gbox = {};
end
if ~isfield(savedata,'eFtilde')
    savedata.eFtilde = [];
end
if ~isfield(savedata,'Pu')
    savedata.Pu = {};
end

% particular solution
if ~isfield(savedata,'G_PU_zero')
    savedata.G_PU_zero = {};
end
if ~isfield(savedata,'G_PU_infty')
    savedata.G_PU_infty = {};
end
if ~isfield(savedata,'PU_A_sum_error')
    savedata.PU_A_sum_error = {};
end

end

% specification-related functions
function doCheck = aux_checkSet(FGunsat,t)
% returns whether the current safe/unsafe set has to be checked for the 
% given time interval, which is not required if the verification has been
% successful in a prior iteration

doCheck = true;

% quick check if set is already fully verified
if isempty(FGunsat)
    doCheck = false;
    return;
end

% check for intersection
if t(2) <= FGunsat(1,1) || t(1) > FGunsat(end,2)
    % upper bound smaller than lower bound of first time interval
    % lower bound larger than upper bound of last time interval
    doCheck = false;
    return;
elseif any(t(1) >= FGunsat(1:end-1,2) & t(2) <= FGunsat(2:end,1))
    % bounds between unverified time intervals
    doCheck = false;
    return;
end

end

function [specUnsat,tFinal] = aux_initSpecUnsat(options,savedata)
% initialize time intervals where all specifications are (not) satisfied
% assume tStart = 0 for now...

tFinal = options.tFinal;
specUnsat = [];

if isfield(savedata,'Funsat') || isfield(savedata,'Gunsat')

    % TODO: shift by start time...
    
    % init whole time horizon with true
    specUnsat = [];
    
    % add unverified time intervals (unsafe sets)
    if isfield(savedata,'Funsat')
        specUnsat = aux_unifySpecUnsat(specUnsat,savedata.Funsat,options.tStart,options.tFinal);
    end
    
    % add unverified time intervals (safe sets)
    if isfield(savedata,'Gunsat')
        specUnsat = aux_unifySpecUnsat(specUnsat,savedata.Gunsat,options.tStart,options.tFinal);
    end
    
    % adjust time horizon if specifications are already satisfied from some
    % time until the end
    tFinal = specUnsat(end,2);
end

end

function specUnsat = aux_unifySpecUnsat(specUnsat,FGunsat,tStart,tFinal)
% specUnsat: double-array nx2
% FGsat:     cell-array mx1 with double-arrays px2

for i=1:length(FGunsat)

    % quick checks
    if isempty(specUnsat)
        % copy FGsat if specUnsat not filled until now
        specUnsat = FGunsat{i};
        continue;
    elseif specUnsat(1,1) == tStart && specUnsat(1,2) == tFinal
        % specUnsat already covers entire time horizon
        break
    end
    
    % loop over all unsat time intervals of current FGsat
    for j=1:size(FGunsat{i},1)
        
        % unsat time interval: add to specUnsat
        t = FGunsat{i}(j,:);
        
        % fully inside
        idx = find(t(1) >= specUnsat(:,1) & t(2) <= specUnsat(:,2),1);
        if ~isempty(idx)
            % already entirely covered -> continue
            continue
        end
        % cases with overlap: cut overlap
        idx = find(t(1) >= specUnsat(:,1) & t(1) <= specUnsat(:,2),1);
        if ~isempty(idx)
            t(1) = specUnsat(idx,2);
        end
        idx = find(t(2) >= specUnsat(:,1) & t(2) <= specUnsat(:,2),1);
        if ~isempty(idx)
            t(2) = specUnsat(idx,1);
        end
        % see if t contains another time interval
        idx = find(t(1) <= specUnsat(:,1) & t(2) >= specUnsat(:,2));
        if ~isempty(idx)
            specUnsat(idx,:) = [];
            % will be implicitly re-added later as t contains it
        end


        % simple cases: no overlap (maximum: point-wise)
        if t(2) < specUnsat(1,1)
            specUnsat = [t; specUnsat];
        elseif t(2) == specUnsat(1,1)
            specUnsat(1,1) = t(1);
        elseif t(1) == specUnsat(end,2)
            specUnsat(end,2) = t(2);
        elseif t(1) > specUnsat(end,2)
            specUnsat = [specUnsat; t];
        else
            idx = find(t(1) >= specUnsat(1:end-1,2) & t(2) <= specUnsat(2:end,1));
            if ~isempty(idx)
                specUnsat = [specUnsat(1:idx,:); t; specUnsat(idx+1:end,:)];
                idxMerge = 0;
                while ~isempty(idxMerge)
                    idxMerge = find(specUnsat(2:end,1) == specUnsat(1:end-1,2));
                    if ~isempty(idxMerge)
                        specUnsat = [specUnsat(1:idx-1,:); [specUnsat(idx,1), specUnsat(idx+1,2)]; specUnsat(idx+2:end,:)];
                    end
                end
            end
        end
        
    end
    
end

end

function [res,set,savedata] = aux_quickCheck(obj,set,timeInterval,...
    Rout_center,eredtotal,Rout_error,options,savedata)

% default value: all ok
res = true;

% check for simple exit (safe set violated)
for i=1:length(savedata.Gunsat)
    
    % propagate generator matrix of GfastInner (for later steps)
    if options.G{i}.fastInner
        C = options.G{i}.set.A;
        % note that G_GfastInner contains mapping by C from i-th safe set,
        % but GfastInner_add does not (so it can be used for all safe sets)
        set.G_GfastInner{i} = set.G_GfastInner{i} + sum(abs(C * set.GfastInner_add),2);
        
        % perform actual check
        if aux_checkSet(savedata.Gunsat{i},timeInterval)
            d = options.G{i}.set.P.b;
            Grem = sum(abs([C * obj.C * generators(set.enc), ...
                C * obj.C * diag(set.boxFstartset_Gbox + set.FtildeuTrans_Gbox + set.G_PUtotal_infty)]),2);
            
            % previously: max(-d + C*c + sum(abs(C*G),2)) > Rout_error(k)
            % but now sum(abs(C*G),2) aims to avoid re-computations
            checkval = max(-d + C*Rout_center + set.G_GfastInner{i} + Grem);
            if checkval > Rout_error - options.errR2Y * eredtotal
                % inner-approximation will be outside of safe set
                res = false; return
            elseif checkval < 0
                % outer-approximation inside safe set -> time interval ok
                savedata.Gunsat{i} = aux_removeFromUnsat(savedata.Gunsat{i},timeInterval);
            end
            
        end
        
    end
    
end

% check for simple exit (unsafe set violated)
for i=1:length(savedata.Funsat)
    
    % propagate generator matrix of GfastInner (for later steps)
    if options.F{i}.fastInner
        C = options.F{i}.set.A;
        % note that G_GfastInner contains mapping by C from i-th safe set,
        % but GfastInner_add does not (so it can be used for all safe sets)
        set.F_GfastInner{i} = set.F_GfastInner{i} + sum(abs(C * set.GfastInner_add),2);
        
        % perform actual check
        if aux_checkSet(savedata.Funsat{i},timeInterval)
            d = options.F{i}.set.P.b;
            Grem = sum(abs([C * obj.C * generators(set.enc), ...
                C * obj.C * diag(set.boxFstartset_Gbox + set.FtildeuTrans_Gbox + set.G_PUtotal_infty)]),2);
            
            % previously: max(-d + C*c + sum(abs(C*G),2)) > Rout_error(k)
            % but now sum(abs(C*G),2) aims to avoid re-computations
            checkval = max(-d + C*Rout_center + set.G_GfastInner{i} + Grem);
            if checkval > Rout_error - options.errR2Y * eredtotal
                % inner-approximation will intersect unsafe set
                res = false; return
            elseif checkval < 0
                % outer-approximation does not intersect unsafe set ->
                % time interval is verified
                savedata.Funsat{i} = aux_removeFromUnsat(savedata.Funsat{i},timeInterval);
            end
            
        end
        
    end
    
end

end

function FGunsat = aux_removeFromUnsat(FGunsat,t)
% adapt FGunsat so that timeInterval is not part of time intervals covered

FGunsat_col = reshape(FGunsat',numel(FGunsat),1);
if mod(sum(t(1) >= FGunsat_col),2) ~= 0
    % lower bound starts inside unverified time interval
    % t(1) \in [ timeInterval )
    idx = find(t(1) >= FGunsat(:,1) & t(1) <= FGunsat(:,2));
    
    if t(2) <= FGunsat(idx,2)
        if t(1) > FGunsat(idx,1)
            FGunsat = [FGunsat(1:idx-1,:); [FGunsat(idx,1), t(1)]; FGunsat(idx:end,:)];
            idx = idx + 1;
        end
        % split, potential merge later
        FGunsat(idx,1) = t(2);
        t = [];
    else
        % remove interval, potential merge later
        FGunsat(idx,2) = t(1);
        if idx < size(FGunsat,1)
            t(1) = FGunsat(idx+1,1);
        end
        if t(2) <= t(1)
            t = [];
        end
    end
end

% now: lower bound starts in between unverified time intervals or at
% maximum at the start point of an unverified set
% t(1) \in [ notTimeInterval )
while ~isempty(t)
    idx = find(t(1) <= FGunsat(:,1),1,'first');
    % upper bound is at least at the beginning of the next time interval
    if t(2) <= FGunsat(idx,2)
        % split, potential merge later
        FGunsat(idx,1) = t(2);
        t = [];
    else
        % remove entire thing (full time interval verified)
        if idx < size(FGunsat,1)
            t(1) = FGunsat(idx,2);
            if t(2) < FGunsat(idx+1,1)
                t = [];
            end
        else
            t = [];
        end
        FGunsat(idx,:) = [];
    end
end

% remove
idxRemove = abs(FGunsat(:,2) - FGunsat(:,1)) < 1e-14;
FGunsat(idxRemove,:) = [];

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
	aux_outputset(obj,isU,set,Rcont_error,options,Rcont_tp_error)
% computes the time-point and time-interval output set as well as the
% contained over-approximation error;
% remark: tp ... end of current step, ti ... start to end of current step
% note: homogeneous and particular solutions are first mapped by the output
% equation and only then added (saving memory allocation when the state
% dimension is much larger than the output dimension)

% additional inclusion of check...?

if isscalar(obj.C) && obj.C == 1 && representsa_(options.V,'origin',eps) && ~any(options.vTrans)
    % y = x ... consequently, errors are also equal
    if isU
        Rout = zonotope([center(set.enc) + set.boxFstartset_center + ...
            set.FtildeuTrans_center, ...
            generators(set.enc), set.G_PUtotal_zero, ...
            diag(set.boxFstartset_Gbox + set.FtildeuTrans_Gbox + set.G_PUtotal_infty)]);
    else
        Rout = zonotope([center(set.enc) + set.boxFstartset_center + ...
            set.FtildeuTrans_center, ...
            generators(set.enc), diag(set.boxFstartset_Gbox + set.FtildeuTrans_Gbox)]);
    end
    Rout_error = Rcont_error;
    
    % time-point solution
    if isU
        Rout_tp = zonotope(set.Hstartp.c, [set.Hstartp.G, ...
            set.G_PUtotal_zero, diag(set.G_PUtotal_infty)]);
    else
        Rout_tp = set.Hstartp;
    end
    Rout_tp_error = Rcont_tp_error;
    
else
    
    if ~isscalar(obj.C) && representsa_(options.V,'origin',eps) && ~any(options.vTrans)
        % y = Cx ... errors are scaled
        
        if isU
            Rout = zonotope([obj.C * (center(set.enc) + set.boxFstartset_center + ...
                set.FtildeuTrans_center), ...
                obj.C * generators(set.enc), obj.C * set.G_PUtotal_zero, ...
                obj.C * diag(set.boxFstartset_Gbox + set.FtildeuTrans_Gbox + set.G_PUtotal_infty)]);
        else
            Rout = zonotope([obj.C * (center(set.enc) + set.boxFstartset_center + ...
                set.FtildeuTrans_center), ...
                obj.C * generators(set.enc), obj.C * diag(set.boxFstartset_Gbox + set.FtildeuTrans_Gbox)]);
        end
        Rout_error = options.errR2Y * Rcont_error;
        
        % time-point solution
        if isU
            Rout_tp = obj.C * zonotope(set.Hstartp.c, [set.Hstartp.G, ...
                set.G_PUtotal_zero, diag(set.G_PUtotal_infty)]);
        else
            Rout_tp = obj.C * set.Hstartp;
        end
        Rout_tp_error = options.errR2Y * Rcont_tp_error;
        

    else
        % general case: y = Cx + V
        
        % compute output set and contained over-approximation error
        if isU
            Rout = zonotope([obj.C * (center(set.enc) + set.boxFstartset_center + ...
                set.FtildeuTrans_center) + options.vTrans, ...
                obj.C * generators(set.enc), obj.C * set.G_PUtotal_zero, ...
                obj.C * diag(set.G_PUtotal_infty + set.boxFstartset_Gbox + set.FtildeuTrans_Gbox), ...
                generators(options.V)]);
        else
            Rout = zonotope([obj.C * (center(set.enc) + set.boxFstartset_center + ...
                set.FtildeuTrans_center) + options.vTrans, ...
                obj.C * generators(set.enc), ...
                obj.C * diag(set.boxFstartset_Gbox + set.FtildeuTrans_Gbox), ...
                generators(options.V)]);
        end
        Rout_error = options.errR2Y * Rcont_error;
        
        % time-point solution
        if isU
            Rout_tp = zonotope([obj.C * center(set.Hstartp) + options.vTrans ...
                obj.C * generators(set.Hstartp), ...
                obj.C * set.G_PUtotal_zero, ...
                obj.C * diag(set.G_PUtotal_infty), generators(options.V)]);
        else
            Rout_tp = obj.C * set.Hstartp + options.V + options.vTrans;
        end
        Rout_tp_error = options.errR2Y * Rcont_tp_error;
        
    end
end


end

% error computation
function errval = aux_errOp(S)
% output of error operator err(*): computes radius of smallest hypersphere
% centered at the origin(!), which encloses S
% this function represents the predominant error measure in this algorithm
% currently not used, because call to costly... instead, the values are
% computed directly in the respective functions

if isa(S,'zonotope')
    % faster computation than calling the interval-constructor
    errval = vecnorm( sum(abs(generators(S)),2) + abs(center(S)) );
elseif isa(S,'interval')
    errval = vecnorm(max(-infimum(S),supremum(S)));
else
    % for all others, convert to interval
    temp = interval(S);
    errval = radius(or(-temp,temp));
end

end

function ebar = aux_ebar(t,k,e,ebar,timeStep,tFinal)
% computes the error bounds for all (accumulating, non-accumulating, and
% reduction) errors for the current step (only the allowed increase in the
% current step! this is different from the sum until the current time)

% compute combined error margin for eacc and enonacc this step: given by
% the maximum error minus the final value for the reduction error and the
% sum of all accumulating errors until the current step
if k == 1
    ebar.remacc(k,1) = e.emax - ebar.red_e(end);
else
    ebar.remacc(k,1) = e.emax - ebar.red_e(end) - e.acctotal(k-1);
end

% the accumulating error bound can be at most the value of the combined
% error margin scaled by the fraction the time step size covers of the
% remaining time
ebar.acc(k,1) = ebar.remacc(k,1) * timeStep / (tFinal - t);


% read out value for admissible reduction error from curve

% find index of nearest time value
t_idx = find(t + timeStep <= ebar.red_t,1,'first');
% use the corresponding linearly interpolated error value as the admissible
% error (also, subtract the reduction error until now, since the
% adaptive-reduce function only considers *local* bounds)
ebarredt_max = ebar.red_e(t_idx-1) + ...
        (ebar.red_e(t_idx) - ebar.red_e(t_idx-1)) * ...
        ( t + timeStep - ebar.red_t(t_idx-1) ) / ...
        ( ebar.red_t(t_idx) - ebar.red_t(t_idx-1) );
if k == 1
    ebar.red(k,1) = ebarredt_max;
    ebar.rem(k,1) = e.emax - ebarredt_max;
else
    ebar.red(k,1) = ebarredt_max - e.redtotal(k-1);
    ebar.rem(k,1) = e.emax - ebarredt_max - e.acctotal(k-1);
end

end

function [e,set,expmat] = aux_errors(fullcomp,obj,isU,isuconst,G_U,u,...
    set,e,expmat,timeStep,timeStepIdx,savedata)
% computes all errors, that is,
%   accumulating:     e.PU_tkplus1,
%   non-accumulating: e.linComb, e.PU_tauk e.F, e.Ftilde,
%   and reduction:    e.red (assigned to 0 as reduction is performed later)
% for the given time step size
% furthermore, auxiliary variables such as F, Ftilde are computed, as well
% as the set PU which is closely related to its respective error

% compute exponential matrix
if timeStepIdx == 0
    % new time step size
    expmat.Deltatk = expm(obj.A * timeStep);
else
    % read from memory
    expmat.Deltatk = savedata.expmatDeltatk{timeStepIdx};
end
expmat.conv = true;

% compute accumulating error
% - consists of eps_PU_tkplus1

% compute PU and eps_PU_tkplus1 (eps_PU_tauk only if fullcomp, see below)
if ~isU
    set.eAtkPU = zeros(obj.dim,1);
    e.PU_tkplus1 = 0;
elseif timeStepIdx == 0
    [set.G_PU_zero,set.G_PU_infty,set.PU_A_sum_error,expmat] = aux_PU(obj,G_U,expmat,timeStep);
    % if sums in computation did not converge, exit now
    if ~expmat.conv
        return;
    end
    % over-approximation of the Hausdorff distance between the exact
    % particular solution PU and the computed one
    e.PU_tkplus1 = vecnorm( sum(abs( (expmat.tk * set.PU_A_sum_error) * G_U ),2) ) ...
        + vecnorm( sum(abs( expmat.tk * set.G_PU_infty ),2) );
    % note: G_PU_zero, G_PU_infty, PU_A_sum_error are saved outside
else
    % re-use values from memory (as time step size already used)
    set.G_PU_infty = savedata.G_PU_infty{timeStepIdx};
    set.PU_A_sum_error = savedata.PU_A_sum_error{timeStepIdx};
    e.PU_tkplus1 = vecnorm( sum(abs( (expmat.tk * set.PU_A_sum_error) * G_U ),2) ) ...
        + vecnorm( sum(abs( expmat.tk * set.G_PU_infty ),2) );
    set.G_PU_zero = timeStep * G_U;
end


% compute non-accumulating error
if fullcomp
% - consists of eps_linComb + 2*err(F H*) + 2*err(Ftilde u) + eps_PU_tauk

% eps_linComb
G_minus = (expmat.Deltatk - eye(obj.dim)) * generators(set.startset);
if strcmp(obj.name,'iss')
    % much faster and similarly precise for ISS system
    % TODO: investigate why this is the case and generalize
    e.linComb = sqrt( size(G_minus,2) ) * (sqrt(norm(G_minus,1) * norm(G_minus,Inf)));
else
    e.linComb = sqrt( size(G_minus,2) ) * norm(G_minus,2);
end
% note: if 2-norm too costly (only ISS for now), switch to:
% 	- norm(G_minus,'fro') ... Frobenius norm > 2-norm
%   - sqrt(norm(G_minus,1)*norm(G_minus,Inf)) ... always > 2-norm
% (... but then error becomes larger so step sizes become smaller)

% curvature errors (time-interval error in state and input)
if timeStepIdx == 0 || ~savedata.fullcomp(timeStepIdx)
    [set,expmat] = aux_F_Ftilde(obj,u,set,expmat,timeStep);
    % if sums in computation did not converge, exit now
    if ~expmat.conv
        return;
    end
    % (F, Fcenter, Frad, Ftilde, Ftildecenter, and Ftilderad saved outside)
else
    set.F = savedata.F{timeStepIdx};
    set.Fcenter = savedata.Fcenter{timeStepIdx};
    set.Frad = savedata.Frad{timeStepIdx};
    set.Ftilde = savedata.Ftilde{timeStepIdx};
    set.Ftildecenter = savedata.Ftildecenter{timeStepIdx};
    set.Ftilderad = savedata.Ftilderad{timeStepIdx};
end
% original computation:
% - e.F = 2 * aux_errOp(set.F * set.startset);
% following version to increase computational efficiency
set.boxFstartset_center = set.Fcenter * set.startset.c;
Fstartset_G = [set.Fcenter * set.startset.G, ...
    diag(set.Frad * sum(abs([set.startset.c, set.startset.G]),2))];
set.boxFstartset_Gbox = sum(abs(Fstartset_G),2);
e.F = 2 * vecnorm(set.boxFstartset_Gbox + abs(set.boxFstartset_center));


% original computation:
% - set.FtildeuTrans = set.Ftilde * zonotope(u);
% - e.Ftilde = 2 * aux_errOp(set.FtildeuTrans);
% following version to increase computational efficiency
set.FtildeuTrans_center = zeros(obj.dim,1);
set.FtildeuTrans_Gbox = zeros(obj.dim,1);
e.Ftilde = 0;
if any(u)
    if isuconst && timeStepIdx > 0 && savedata.fullcomp(timeStepIdx)
        set.FtildeuTrans_center = savedata.FtildeuTrans_center{timeStepIdx};
        set.FtildeuTrans_Gbox = savedata.FtildeuTrans_Gbox{timeStepIdx};
        e.Ftilde = savedata.eFtilde{timeStepIdx};
    else
        if timeStepIdx > 0 && savedata.fullcomp(timeStepIdx)
            % read Ftildecenter and Ftilderad from memory, but compute
            % Ftilde * uTrans new (since u might have changed)
            set.Ftildecenter = savedata.Ftildecenter{timeStepIdx};
            set.Ftilderad = savedata.Ftilderad{timeStepIdx};
        end
        set.FtildeuTrans_center = set.Ftildecenter * u;
        set.FtildeuTrans_Gbox = set.Ftilderad * abs(u);
        % quicker version (removing unnecessary functions w.r.t errOp-call)
        e.Ftilde = 2 * vecnorm(set.FtildeuTrans_Gbox + abs(set.FtildeuTrans_center));
    end
end


% err(e^At PU)
% over-approximation of the Hausdorff distance between the computed
% particular solution PU and 0 (minimum solution over time interval)
e.PU_tauk = 0;
if isU
    e.PU_tauk = vecnorm( sum(abs(expmat.tk * set.G_PU_zero),2) + ...
        sum(abs(expmat.tk * set.G_PU_infty),2) );
    % ... equal to aux_errOp(set.eAtkPU)
end

% total non-accumulating error
e.nonacc = e.linComb + e.PU_tauk + e.F + e.Ftilde;

else

% use NaN values for completeness
e.F = NaN;
e.Ftilde = NaN;
e.linComb = NaN;
e.PU_tauk = NaN;
e.nonacc = NaN;

end


% total accumulating error (contribution of current step)
% note that e.PU_tauk over-approximates e.PU_tkplus1, but for the
% accumulation to the next step, e.PU_tkplus1 is used
e.acc = e.PU_tkplus1;

% init value for reduction error
e.red = 0;

end

function [set,expmat] = aux_F_Ftilde(obj,u,set,expmat,timeStep)
% computation of F and Ftilde for a given time step size; currently, these
% variables are computed by a Taylor series until floating-point precision,
% i.e., we increase the truncation order until the additional values are so
% small that the stored number (finite precision!) does not change anymore

% skip computation of Ftilde if u is all-zero

% load data from object/options structure
A = obj.A;
n = obj.dim;

% is there an input vector?
isu = any(u);

% initialize auxiliary variables and flags for loop
Asum_pos_F = zeros(n);
Asum_neg_F = zeros(n);
stoploop_F = false;
if isu
    Asum_pos_Ftilde = zeros(n);
    Asum_neg_Ftilde = zeros(n);
    stoploop_Ftilde = false;
else
    set.Ftilde = [];
    set.Ftildecenter = [];
    set.Ftilderad = [];
    stoploop_Ftilde = true;
end

eta = 1;
while true
    % exponential: 1:eta
    
    % compute powers
    expmat.Apower = aux_getApower(eta,A,expmat.Apower);
    Apower_eta = expmat.Apower{eta};
    
    % F starts at eta = 2, so skip for eta = 1
    if eta==1; eta = eta + 1; continue; end
    
    % tie/inputTie: 2:eta
    % note: usually, inputTie goes to eta+1 (with eta from F), but since we
    % compute terms until floating-point precision, this does not need to
    % be respected (only if we were to use a remainder term E, which then
    % would necessarily need to be adapted to a fixed eta)
    
    % compute factor (factorial already included in powers of A)
    exp1 = -(eta)/(eta-1); exp2 = -1/(eta-1);
    factor = ((eta)^exp1-(eta)^exp2) * timeStep^eta; % previously: /factorial(eta)
    
    if ~stoploop_F
        [expmat.Apos,expmat.Aneg] = ...
            aux_getAposneg(eta,n,expmat.Apos,expmat.Aneg,Apower_eta);
        
        % if new term does not change result anymore, loop to be finished
        Asum_add_pos_F = factor*expmat.Aneg{eta};
        Asum_add_neg_F = factor*expmat.Apos{eta};
        
        % safety check (if time step size too large, then the sum converges
        % too late so we already have Inf values)
        % previous criterion: any(any(isinf(Asum_add_pos_F))) || any(any(isinf(Asum_add_neg_F)))
        % ... but very costly for large A!
        if eta == 75
            set.F = []; set.Ftilde = [];
            expmat.conv = false; return
%             throw(MException('reach_adaptive:notconverging',...
%                 'Time Step Size too big for computation of F.'));
        end
        
        % compute ratio for floating-point precision
        if all(all(Asum_add_pos_F <= eps * Asum_pos_F)) && ...
                all(all(Asum_add_neg_F >= eps * Asum_neg_F))
            stoploop_F = true;
        end

        % compute powers; factor is always negative
        Asum_pos_F = Asum_pos_F + Asum_add_pos_F; 
        Asum_neg_F = Asum_neg_F + Asum_add_neg_F;

    end
    
    if ~stoploop_Ftilde
        [expmat.Apos,expmat.Aneg] = ...
            aux_getAposneg(eta-1,n,expmat.Apos,expmat.Aneg,expmat.Apower{eta-1});
        
        % if new term does not change result anymore, loop to be finished
        % we require one additional division by eta as the terms in expmat
        % are divided by (eta-1)! instead of eta! as required
        Asum_add_pos_Ftilde = factor*expmat.Aneg{eta-1} / eta;
        Asum_add_neg_Ftilde = factor*expmat.Apos{eta-1} / eta;
        
        % safety check (if time step size too large, then the sum converges
        % too late so we already have Inf values)
        % previous criterion: any(any(isinf(Asum_add_pos_Ftilde))) || any(any(isinf(Asum_add_neg_Ftilde)))
        % ... but very costly for large A!
        if eta == 75
            set.F = []; set.Ftilde = [];
            expmat.conv = false; return
%             throw(MException('reach_adaptive:notconverging',...
%                 'Time Step Size too big for computation of Ftilde.'));
        end
        
        % compute ratio for floating-point precision
        if all(all(Asum_add_pos_Ftilde <= eps * Asum_pos_Ftilde)) && ...
                all(all(Asum_add_neg_Ftilde >= eps * Asum_neg_Ftilde)) 
            stoploop_Ftilde = true;
        end

        % compute powers; factor is always negative
        Asum_pos_Ftilde = Asum_pos_Ftilde + Asum_add_pos_Ftilde; 
        Asum_neg_Ftilde = Asum_neg_Ftilde + Asum_add_neg_Ftilde;
    end
    
    % instantiate interval matrices if converged
    if stoploop_F
        set.F = interval(Asum_neg_F,Asum_pos_F);
        set.Frad = 0.5*(Asum_pos_F - Asum_neg_F);
        set.Fcenter = Asum_neg_F + set.Frad;
    end
    if stoploop_Ftilde && isu
        set.Ftilde = interval(Asum_neg_Ftilde,Asum_pos_Ftilde);
        set.Ftilderad = 0.5*(Asum_pos_Ftilde - Asum_neg_Ftilde);
        set.Ftildecenter = Asum_neg_Ftilde + set.Ftilderad;
    end

    % exit loop if both converged
    if stoploop_F && stoploop_Ftilde
        break;
    end
    
    % increment eta
    eta = eta + 1;
end

end

function [Apower,Apower_abs] = aux_getApower(eta,A,Apower,A_abs,Apower_abs)
% this function ensures that the eta-th power of A and |A| is computed
% (this is necessary, since we do not know the largest power in advance,
% and we want to save computations as much as possible)
% we do not compute A^eta but A^eta / eta! instead to increase the stability
% -> this has to be taken into account in all use cases!
% this is currently not enacted for |A|^eta


% check A^eta
if length(Apower) >= eta
    % read from memory
else
    % compute all terms A^i/i! until eta
    maxeta = length(Apower);
    for i=maxeta:eta-1
        Apower{i+1} = Apower{i} * A / (i+1);
        % sparse/full  storage for more efficiency
        if nnz(Apower{i+1}) / (size(Apower{i+1},1)^2) < 0.1
            Apower{i+1} = sparse(Apower{i+1});
        else
            Apower{i+1} = full(Apower{i+1});
        end
    end
end

if nargout == 2
    % check |A|^eta
    if length(Apower_abs) >= eta
        % read from memory
    else
        % compute all powers |A|^i until eta
        maxeta = length(Apower_abs);
        for i=maxeta:eta-1
            Apower_abs{i+1} = Apower_abs{i}*A_abs;
            % sparse/full storage for more efficiency
            if nnz(Apower_abs{i+1}) / (size(Apower_abs{i+1},1)^2) < 0.1
                Apower_abs{i+1} = sparse(Apower_abs{i+1});
            else
                Apower_abs{i+1} = full(Apower_abs{i+1});
            end
        end
    end
end

end

function [Apos,Aneg] = aux_getAposneg(eta,n,Apos,Aneg,Apower_eta)
% the separation of A^eta into positive and negative indices can be
% precomputed and saved; Apower_eta has to match eta correctly!

if length(Apos) >= eta && ~isempty(Apos{eta})
    % ... then also length(Aneg) >= eta
    % read from memory
    
else
    
    % old method (same as in tie)
%     % obtain positive and negative parts
%     pos_ind = Apower_eta > 0;
%     neg_ind = Apower_eta < 0;
% 
%     % init positive and negative parts at index eta
%     Apos{eta} = zeros(n);
%     Aneg{eta} = zeros(n);
%     
%     % insert values from A^eta
%     Apos{eta}(pos_ind) = Apower_eta(pos_ind);
%     Aneg{eta}(neg_ind) = Apower_eta(neg_ind);
    
    % new method (slightly faster)
    Aneg{eta} = Apower_eta;
    Apos{eta} = Apower_eta;
    Aneg{eta}(Aneg{eta} > 0) = 0;
    Apos{eta}(Apos{eta} < 0) = 0;
end

end

% particular solutions
function [Pu,expmat] = aux_Pu(obj,u,expmat,timeStep)
% computation of the particular solution due to the input vector u using
% a Taylor series where the truncation order is increased until the
% additional values are so small that the stored number (finite precision!)
% does not change anymore; in case the inverse of A exists, we directly
% compute the analytical solution (where the exponential matrix is also
% only computed until finite precision)

if ~any(u)
    Pu = zeros(obj.dim,1);

elseif expmat.isAinv
    Pu = expmat.Ainv * (expmat.Deltatk - eye(obj.dim)) * u;
    
else    
    % compute by sum until floating-point precision (same as for PU)
    % formula: \sum_{j=1}^\infty \frac{A^{j-1}}{j!} timeStep^{j}
    
    % initialize truncation order
    eta = 1;
    
    % first term
    Asum = timeStep * eye(obj.dim);
    
    % loop until Asum no longer changes (additional values too small)
    while true
        % increment truncation order
        eta = eta + 1;
        
        % get A^eta-1
        expmat.Apower = aux_getApower(eta-1,obj.A,expmat.Apower);
        Apower_etaminus1 = expmat.Apower{eta-1};
        
        % compute additional term (division by (eta-1)! already included in
        % Apower_etaminus1, so one additional /eta required)
        addTerm = Apower_etaminus1 / eta * timeStep^eta;
        
        % safety check (if time step size too large, then the sum converges
        % too late so we already have Inf values)
        if any(any(isinf(addTerm)))
            % this error is currently not handled in calling function
            throw(MException('reach_adaptive:notconverging',...
                'Time Step Size too big for computation of Pu.'));
        end
        
        % if new term does not change stored values in Asum, i.e., all
        % entries are below floating-point accuracy -> stop loop
        if all(all(abs(addTerm) <= eps * abs(Asum)))
            break;
        end
        
        % add term to current Asum
        Asum = Asum + addTerm;
    end
    
    % compute particular solution due to input vector
    Pu = Asum * u;
end

end

function [G_PU_zero,G_PU_infty,A_sum_error,expmat] = aux_PU(obj,G_U,expmat,timeStep)
% computation of the particular solution due to the uncertain input set U
% as well as the contained over-approximation error with respect to the
% exact particular solution; using a Taylor series where the truncation
% order is increased until the additional values are so small that the
% stored number (finite precision!) does not change anymore
% we use only the generator matrix to spare some call to zonotope functions

% initialize particular solution (only generator matrices)
G_PU_zero = timeStep * G_U;
PU_diag = sum(abs(G_PU_zero),2);

% initialize errors
G_PU_infty = [];
A_sum_error = zeros(obj.dim);

A = obj.A;

% loop until floating-point precision
stoploop = false;
eta = 1;
while true
    
    % compute powers of A
    expmat.Apower = aux_getApower(eta,A,expmat.Apower);
    Apower_eta = expmat.Apower{eta};
    
    % additional term (Apower_eta already contains division by (eta)!, thus
    % we require one more /(eta+1) to get correct denominator)
    addG_PU = Apower_eta / (eta+1) * timeStep^(eta+1) * G_U;
    addG_PU_diag = sum(abs(addG_PU),2);
    
    % safety check (if time step size too large, then the sum converges
    % too late so we already have Inf values)
    % previous criterion: any(any(isinf(addPU_diag))) || any(any(isnan(addPU_diag)))
    % ... but very costly for large A!
    if eta == 75
        expmat.conv = false; return
%         throw(MException('reach_adaptive:notconverging',...
%             'Time Step Size too big for computation of PU.'));
    end
    
    % check if floating-point precision reached
    if all( abs(addG_PU_diag) <= eps * abs(PU_diag) )
        stoploop = true;
    end
    
    % add term to simplified value for convergence
    PU_diag = PU_diag + addG_PU_diag;
    
    % error terms
    A_sum_error = A_sum_error + Apower_eta / (eta+1) * timeStep^(eta+1);
    G_PU_infty = [G_PU_infty, addG_PU];
    
    % break loop
    if stoploop
%         disp("eta: " + eta);
        break;
    end
    
    % increment eta
    eta = eta + 1;

end

end

function [set,e] = aux_reduce(obj,k,set,e,ebar,expmat)
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

% state dimension
n = size(expmat.Deltatk,1);

% the _infty parts are always boxed (which induces no additional error),
% saved as a vector... correct handling in computation of output set
G_eAtkPU_infty = sum(abs( expmat.tk * set.G_PU_infty ),2);
set.G_PUtotal_infty = set.G_PUtotal_infty + G_eAtkPU_infty;

% generator matrix of additional solution e^Atk PU_zero
G_eAtkPU_zero = expmat.tk * set.G_PU_zero;
Gabs_eAtkPU_zero = abs(G_eAtkPU_zero);

% for fast falsification
set.GfastInner_add = obj.C * G_eAtkPU_zero;

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

if isempty(omega_max_min) || omega_max_min > ebar.red(k)
    % no generators can be reduced
    e.red = 0;
    set.nrG_PUtotal_zero = set.nrG_PUtotal_zero + nrG_eAtkPU_zero;
else
    set.nrG_PUtotal_zero = size(set.girard_PUtotal_zero,2);
    [~,idx] = mink(set.girard_PUtotal_zero,set.nrG_PUtotal_zero);
    
    % compute over-approximation of dH (omega_max)
    sum_temp = zeros(n,1);
    e.red = 0; redIdx = set.nrG_PUtotal_zero;
    for ijk=1:set.nrG_PUtotal_zero
        sum_temp = sum_temp + set.Ghat_PUtotal_zero(:,idx(ijk));
        omega_max_ijk = 2 * vecnorm(sum_temp,2);
        if omega_max_ijk > ebar.red(k)
            redIdx = ijk-1; break;
        else
            e.red = omega_max_ijk;
        end
    end
    
%     % version using cumsum (usually slower as redidx usually close to 1)
%     auxval = cumsum(set.Ghat_PUtotal_zero(:,idx),2);
%     omega_max = 2 * vecnorm(auxval,2); % faster than: sqrt(sum(auxval.^2,1))
%     % index until which gens are reduced
%     redIdx = find(omega_max <= ebar.red(k),1,'last');
%     % maximum reduction error which can be induced in this step
%     e.red = omega_max(redIdx);
    
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

function [set,e] = aux_reduce_rad(obj,k,set,e,ebar,expmat)
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
% --- saved in set-struct for usage in successive steps
% ... G_PUtotal_zero nevers contains any girard-0-generators

% system dimension
n = size(expmat.Deltatk,1);

% the _infty parts are always boxed (which induces no additional error),
% saved as a vector... correct handling in computation of output set
G_eAtkPU_infty = sum(abs( expmat.tk * set.G_PU_infty ),2);
set.G_PUtotal_infty = set.G_PUtotal_infty + G_eAtkPU_infty;

% generator matrix of additional solution e^Atk PU_zero
G_eAtkPU_zero = expmat.tk * set.G_PU_zero;
Gabs_eAtkPU_zero = abs(G_eAtkPU_zero);

% for fast falsification
set.GfastInner_add = obj.C * G_eAtkPU_zero;

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

% reconstruct matrices for PUtotal_zero = PU_zero(tk) + e^Atk PU_zero(Delta tk)
% pre-allocation version
% if set.nrG_PUtotal_zero+nrG_eAtkPU_zero > set.nrG_PUtotal_zero_full
%     set.nrG_PUtotal_zero_full = set.nrG_PUtotal_zero_full + set.nrG_PUtotal_zero_incr;
%     set.G_PUtotal_zero(:,end+1:end+set.nrG_PUtotal_zero_incr) = zeros(n,set.nrG_PUtotal_zero_incr);
%     set.Gabs_PUtotal_zero(:,end+1:end+set.nrG_PUtotal_zero_incr) = zeros(n,set.nrG_PUtotal_zero_incr);
% end
set.G_PUtotal_zero(:,set.nrG_PUtotal_zero+1:set.nrG_PUtotal_zero+nrG_eAtkPU_zero) = G_eAtkPU_zero;
set.Gabs_PUtotal_zero(:,set.nrG_PUtotal_zero+1:set.nrG_PUtotal_zero+nrG_eAtkPU_zero) = Gabs_eAtkPU_zero;

% reconstruct full girard metric
set.girard_PUtotal_zero(:,set.nrG_PUtotal_zero+1:set.nrG_PUtotal_zero+nrG_eAtkPU_zero) = girard_eAtkPU_zero(~axisAligned);

% compute over-approximation of dH (omega_rad): only first value
% check whether this value is already larger than reduction error bound
[~,minidx] = min(set.girard_PUtotal_zero);
omega_rad_min = 2 * vecnorm(set.Gabs_PUtotal_zero(:,minidx),2);

if isempty(omega_rad_min) || omega_rad_min > ebar.red(k)
    % no generators can be reduced
    e.red = 0;
    set.nrG_PUtotal_zero = set.nrG_PUtotal_zero + nrG_eAtkPU_zero;
else
    % linear search to find index
    set.nrG_PUtotal_zero = size(set.girard_PUtotal_zero,2);
    [~,idx] = mink(set.girard_PUtotal_zero,set.nrG_PUtotal_zero);
    
    % compute over-approximation of dH (omega_rad)
    sum_temp = zeros(n,1);
    e.red = 0;
    for ijk=1:set.nrG_PUtotal_zero
        sum_temp = sum_temp + set.Gabs_PUtotal_zero(:,idx(ijk));
        omega_rad_ijk = vecnorm(sum_temp,2);
        if omega_rad_ijk > ebar.red(k)
            redIdx = ijk-1; break;
        else
            e.red = omega_rad_ijk;
        end
    end
    
    % version using cumsum (usually slower as redidx usually close to 1)
%     omega_rad = vecnorm(cumsum(set.Gabs_PUtotal_zero(:,idx),2),2);
    % index until which gens are reduced
%     redIdx = find(omega_rad <= ebar.red(k),1,'last');
    % maximum reduction error which can be induced in this step
%     e.red = omega_rad(redIdx);
    
    % indices which are reduced
    idxRed = idx(1:redIdx);

    % add box made of generators selected for reduction to infty part
    set.G_PUtotal_infty = set.G_PUtotal_infty + sum(set.Gabs_PUtotal_zero(:,idxRed),2);

    % remove indices which were selected for reduction
    set.girard_PUtotal_zero(:,idxRed) = [];
    set.Gabs_PUtotal_zero(:,idxRed) = [];
    set.G_PUtotal_zero(:,idxRed) = [];
    set.nrG_PUtotal_zero = set.nrG_PUtotal_zero - redIdx;
%     set.nrG_PUtotal_zero_full = set.nrG_PUtotal_zero_full - redIdx;
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

function coeff = aux_coeff(cnti,fullcomp,coeff,e,timeStep,cnt)
% compute the coefficients for the approximation functions which model the
% behavior of the individual error over Delta t

% note the order of the individual terms over Delta t
% - eps_linComb:    linear
% - eps_F:          quadratic
% - eps_Ftilde:     quadratic
% - eps_PU:         quadratic

% coefficients for approximation functions
% - linear:         e(Delta t) =                 b * Delta t
% - quadratic:      e(Delta t) = a * Delta t^2 + b * Delta t
% ... no offset since all errors go to 0 for delta t -> 0

if cnti == 1
    % initialize coefficients (guessing that b = 0)
    % TODO: check if re-use from last step possible
    
    if fullcomp
    % linComb and time-interval error of PU are linear approximation
    % functions
    coeff.linComb.b = e.linComb / timeStep;
    coeff.PU_tauk.b = e.PU_tauk / timeStep;

    % all others are quadratic approximation functions
    % here, we guess that b = 0
    coeff.F.a = e.F / timeStep^2;
    coeff.F.b = 0;
    coeff.Ftilde.a = e.Ftilde / timeStep^2;
    coeff.Ftilde.b = 0;
    end
    
    coeff.PU_tkplus1.a = e.PU_tkplus1 / timeStep^2;
    coeff.PU_tkplus1.b = 0;
    
else
    % update coefficients
    
    if fullcomp
    % linear approximation function -> take newest value
    coeff.linComb.b = e.linComb / timeStep;
    coeff.PU_tauk.b = e.PU_tauk / timeStep;
    end
    
    % quadratic approximation functions -> take most recent two values
    % a * lastTimeStep^2 + b * lastTimeStep = lastEps
    % a * timeStep^2     + b * timeStep     = eps
    % -> Deltatmat * [a;b] = eps ... solve using '\'-operator
    
    % matrix of time step sizes (quadratic and linear)
    Deltatmat = [cnt.timeStep(cnti-1)^2     cnt.timeStep(cnti-1); ...
                 cnt.timeStep(cnti)^2       cnt.timeStep(cnti)];

    % sanity check of singularity of Deltatmat
    if abs(1/cond(Deltatmat)) < eps
        throw(CORAerror('CORA:notConverged','Estimation of time step size'));
    end

    % update coefficient of approximation function for eps_PU
    tmp = Deltatmat \ [cnt.ePU_tkplus1(cnti-1); e.PU_tkplus1];
    coeff.PU_tkplus1.a = tmp(1); coeff.PU_tkplus1.b = tmp(2);
    
    if fullcomp
    % update coefficient of approximation function for eps_F
    tmp = Deltatmat \ [cnt.eF(cnti-1); e.F];
    coeff.F.a = tmp(1); coeff.F.b = tmp(2);
    
    if coeff.F.b < 0
        % then a certain region is smaller than 0 which cannot ever happen
        % -> use only current value and set coeff.F.b to 0
        coeff.F.b = 0;
        coeff.F.a = e.F / timeStep^2;
    end

    % update coefficient of approximation function for eps_Ftilde
    tmp = Deltatmat \ [cnt.eFtilde(cnti-1); e.Ftilde];
    coeff.Ftilde.a = tmp(1); coeff.Ftilde.b = tmp(2);
    end
end

end

function timeStep = aux_predTimeStep_approxfuncs(t,k,fullcomp,coeff,...
    e,ebar,timeStep,maxTimeStep,tFinal)
% predict a time step size which satisfies the error bound ebar.rem
% (comprising the remaining error for e.acc and e.nonacc) using the
% approximation functions for the individual errors and respecting the
% maximum admissible time step size due possible changes in the input
% vector (options.uTransVec)

% timeStep necessary to satisfy eAcc_bound by equating the linearly
% increasing bound and the approximation function

% in the first step, there is a chance that the initial guess is too large
% thus, we use another condition to decrease the initial guess until a
% reasonable value can be found
% (note: this works also when e.nonacc = [])
if k == 1 && any([e.nonacc / ebar.rem, e.acc / ebar.acc] > 1e3)
    timeStep = 0.01 * timeStep;
    return;
end

% safety factor so that we are more likely to satisfy the bounds in case
% the approximation functions underestimate the errors
safetyFactor = 0.90;

% 1. condition: accumulating error needs to satisfy linearly increasing
% bound until the time horizon
% note: if no inputs, then this value becomes Inf (see min-operation below)
timeStep_pred_linBound = safetyFactor * ...
    (1/coeff.PU_tkplus1.a * (ebar.remacc(k) / (tFinal - t) - coeff.PU_tkplus1.b));

timeStep_pred_erem = Inf;
if fullcomp
% 2. condition: both errors need to satisfy erem for the current step
% ...this yields a quadratic equation
temp_a = coeff.PU_tkplus1.a + coeff.F.a + coeff.Ftilde.a;
temp_b = coeff.PU_tkplus1.b + coeff.F.b + coeff.Ftilde.b ...
            + coeff.linComb.b + coeff.PU_tauk.b;
temp_c = -ebar.rem(k);

% compute two solutions of quadratic equation (one of them is < 0)
% choose maximum of predicted timeSteps -> positive solution
timeStep_pred_erem = safetyFactor * max(...
    [1/(2*temp_a) * (-temp_b + sqrt(temp_b^2 - 4*temp_a*temp_c)); ...
     1/(2*temp_a) * (-temp_b - sqrt(temp_b^2 - 4*temp_a*temp_c))] );
end

% update timeStep by minimum of predicted timeSteps so that both error
% bounds are likely to be satisfied; also, the maximum admissible for the
% current step needs to be respected as the input vector has to be constant
% over any single step
timeStep = min([maxTimeStep;timeStep_pred_linBound;timeStep_pred_erem]);

end

function bisect = aux_bisect(k,cnti,bisect,isU,e,ebar,timeStep,cnt)
% update values stored in bisect (used for bisection algorithm to determine
% a suitable time step size satisfying the error bounds)
% errcheck: true (errors have been fulfilled) / false (opposite)

% we always save the error values for a time step size that is too small
% (which is zero by default) and a time step size that is also too small or
% large enough so that the intersection is between the lower and upper bound

if cnti == 1
    % initialize lower bound with 0
    bisect.lb.timeStep_acc = 0;
    bisect.lb.acc = 0;
    bisect.lb.accok = true;
    bisect.lb.acc_perc = 0;
    
    bisect.lb.timeStep_nonacc = 0;
    bisect.lb.nonacc = 0;
    bisect.lb.nonaccok = true;
    bisect.lb.nonacc_perc = 0;
    
    % upper bound is time step size
    bisect.ub.timeStep_acc = timeStep;
    bisect.ub.acc = e.PU_tkplus1;
    bisect.ub.accok = false;
    if cnt.accok(cnti)
        bisect.ub.accok = true;
    end
    bisect.ub.acc_perc = e.acc / ebar.acc(k);
    
    bisect.ub.timeStep_nonacc = timeStep;
    bisect.ub.nonacc = e.nonacc;
    bisect.ub.nonaccok = false;
    if cnt.nonaccok(cnti)
        bisect.ub.nonaccok = true;
    end
    bisect.ub.nonacc_perc = e.nonacc / ebar.rem(k);
    
else

    % determine whether new time step size is new lb or new ub
    if isU
    if ~cnt.accok(cnti) || timeStep > bisect.ub.timeStep_acc
        % current time step size is always new ub if errcheck not ok
        
        if timeStep > bisect.ub.timeStep_acc
            % current ub becomes lb
            bisect.lb.timeStep_acc = bisect.ub.timeStep_acc;
            bisect.lb.acc = bisect.ub.acc;
            bisect.lb.accok = bisect.ub.accok;
            bisect.lb.acc_perc = bisect.ub.acc_perc;
        end
        % assign new ub
        bisect.ub.timeStep_acc = timeStep;
        bisect.ub.acc = e.PU_tkplus1;
        bisect.ub.accok = cnt.accok(cnti);
        bisect.ub.acc_perc = e.acc / ebar.acc(k);
    else % smaller than before and errcheck ok -> new lb
        bisect.lb.timeStep_acc = timeStep;
        bisect.lb.acc = e.PU_tkplus1;
        bisect.lb.accok = cnt.accok(cnti);
        bisect.lb.acc_perc = e.acc / ebar.acc(k);
    end
    end
    
    if ~cnt.nonaccok(cnti) || timeStep > bisect.ub.timeStep_nonacc
        % current time step size is always new ub if errcheck not ok
        
        if timeStep > bisect.ub.timeStep_nonacc
            % current ub becomes lb
            bisect.lb.timeStep_nonacc = bisect.ub.timeStep_nonacc;
            bisect.lb.nonacc = bisect.ub.nonacc;
            bisect.lb.nonaccok = bisect.ub.nonaccok;
            bisect.lb.nonacc_perc = bisect.ub.nonacc_perc;
        end
        bisect.ub.timeStep_nonacc = timeStep;
        bisect.ub.nonacc = e.nonacc;
        bisect.ub.nonaccok = cnt.nonaccok(cnti);
        bisect.ub.nonacc_perc = e.nonacc / ebar.rem(k);
    else % smaller than before and errcheck ok -> new lb
        bisect.lb.timeStep_nonacc = timeStep;
        bisect.lb.nonacc = e.nonacc;
        bisect.lb.nonaccok = cnt.nonaccok(cnti);
        bisect.lb.nonacc_perc = e.nonacc / ebar.rem(k);
    end
end

end

function timeStep = aux_predTimeStep_bisect(t,k,fullcomp,bisect,isU,...
    e,ebar,timeStep,maxTimeStep,tFinal)
% predict a time step size which satisfies the error bound ebar.rem
% (comprising the remaining error for e.acc and e.nonacc) using a weighted
% bisection approach between the current and the previous value (or 0),
% the maximum admissible time step size due possible changes in the input
% vector (options.uTransVec) is also respected

% special handling for first step
if k == 1
    eacctotal = 0;
else
    eacctotal = e.acctotal(k-1);
end

% BISECTION METHOD
% 1. acc errors
timeStep_pred_eacc = Inf;
slope_acc = 0;
if isU
if bisect.ub.accok
    % extrapolate
    ebaraccend_tFinal = (e.emax - ebar.red_e(end)) / tFinal;
    slope_acc = (bisect.ub.acc - bisect.lb.acc) / ...
        (bisect.ub.timeStep_acc - bisect.lb.timeStep_acc);
    
    if slope_acc < ebaraccend_tFinal
        % slope of error bound larger than slope of error (this occurs due
        % to assumed linearity... actually, the slope of the error will
        % increase for larger values due to higher-order terms)
        timeStep_pred_eacc = Inf;
    else
        timeStep_add = ( ebaraccend_tFinal * (t + bisect.lb.timeStep_acc) ...
            - eacctotal - bisect.lb.acc ) / (slope_acc - ebaraccend_tFinal);
        timeStep_pred_eacc = bisect.lb.timeStep_acc + timeStep_add;
    end
else
    % bisection
    if k == 1 % constant in first step as differences very large
        factor = 0.5;
    else % weighted, using information that decrease at least quadratically
        % ...estimate so that bound fulfilled at 95%
        factor = sqrt( (0.95 - bisect.lb.acc_perc) / ...
            (bisect.ub.acc_perc - bisect.lb.acc_perc) );
    end
    timeStep_pred_eacc = (bisect.lb.timeStep_acc + ...
        factor * (bisect.ub.timeStep_acc - bisect.lb.timeStep_acc));
end
end

% 2. nonacc errors
timeStep_pred_enonacc = Inf;
if fullcomp
if bisect.ub.nonaccok
    % extrapolate
    slope_nonacc = (bisect.ub.nonacc - bisect.lb.nonacc) / ...
    	(bisect.ub.timeStep_nonacc - bisect.lb.timeStep_nonacc);
    if isU
        timeStep_add = (e.emax - ebar.red_e(end)*t/tFinal - eacctotal ...
            - bisect.lb.acc - bisect.lb.nonacc ) / ...
            ( slope_acc + slope_nonacc - ebar.red_e(end)/tFinal );
    else 
        % shorter way
        timeStep_add = (e.emax - bisect.lb.nonacc) / slope_nonacc;
    end
    timeStep_pred_enonacc = bisect.lb.timeStep_nonacc + timeStep_add;
else
    % bisection
    if k == 1 % constant in first step as differences very large
        factor = 0.5;
    else % weighted, estimate so that bound fulfilled at 95%
        factor = (0.95 - bisect.lb.nonacc_perc) / ...
            (bisect.ub.nonacc_perc - bisect.lb.nonacc_perc);
    end
    timeStep_pred_enonacc = (bisect.lb.timeStep_nonacc + ...
            factor * (bisect.ub.timeStep_nonacc - bisect.lb.timeStep_nonacc));
end
% for safety... (if slopes misbehave)
if timeStep_pred_enonacc < 0
    timeStep_pred_enonacc = Inf;
end
end

% find minimum
timeStep = min([timeStep_pred_enonacc,timeStep_pred_eacc,maxTimeStep]);

end

function idx = aux_checkTimeStep(timeStep,savedata)
% checks whether auxiliary sets have already been computed for this time
% step size (allows to save a lot of recomputations)

idx = find(abs(savedata.timeStep - timeStep) < 1e-9,1);
if isempty(idx)
    idx = 0; % set 0 since this can be assigned to a list and is not an index
end

end

function [timeStep,timeStepIdx] = aux_adjustTimeStep(k,timeStep,isU,maxTimeStep,savedata)
% adjusts the predicted time step size in order to avoid recomputations
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
    if ~( all(timeStep > savedata.timeStep) || all(timeStep < savedata.timeStep) )
        timeStep_ub = min(savedata.timeStep(savedata.timeStep > timeStep));
        timeStep_lb = max(savedata.timeStep(savedata.timeStep < timeStep));
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
function [compR, maxTimeStepSpec] = aux_fullcomp(t,specUnsat)
% checks whether at the current point in time the full reachable set has to
% be computed since the specification need to be checked; additionally, we
% return the time until the next change in the boolean value of specSAT

% default values
maxTimeStepSpec = Inf;
compR = true;

% no specifications given
if isempty(specUnsat)
    return
end

% check if current time within bounds (but not at the end)
temp = reshape(specUnsat',numel(specUnsat),1);
idx = find(t < temp,1,'first');
% time until next switch
maxTimeStepSpec = temp(idx) - t;
% compute full reachable set?
if mod(idx,2) == 1
	compR = false;
end

end

function [timeStepIdx,savedata] = aux_savedata(fullcomp,isU,isuconst,...
    set,e,expmat,timeStep,savedata)
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
% propagation matrix
savedata.expmatDeltatk{timeStepIdx,1} = expmat.Deltatk;
% affine dynamics
savedata.fullcomp(timeStepIdx,1) = fullcomp;
if fullcomp
    savedata.F{timeStepIdx,1} = set.F;
    savedata.Fcenter{timeStepIdx,1} = set.Fcenter;
    savedata.Frad{timeStepIdx,1} = set.Frad;
    savedata.Ftilde{timeStepIdx,1} = set.Ftilde;
    savedata.Ftildecenter{timeStepIdx,1} = set.Ftildecenter;
    savedata.Ftilderad{timeStepIdx,1} = set.Ftilderad;
    if isuconst
        savedata.FtildeuTrans_center{timeStepIdx,1} = set.FtildeuTrans_center;
        savedata.FtildeuTrans_Gbox{timeStepIdx,1} = set.FtildeuTrans_Gbox;
        savedata.eFtilde{timeStepIdx,1} = e.Ftilde;
    end
else
    savedata.F{timeStepIdx,1} = [];
    savedata.Fcenter{timeStepIdx,1} = [];
    savedata.Frad{timeStepIdx,1} = [];
    savedata.Ftilde{timeStepIdx,1} = [];
    savedata.Ftildecenter{timeStepIdx,1} = [];
    savedata.Ftilderad{timeStepIdx,1} = [];
    if isuconst
        savedata.FtildeuTrans_center{timeStepIdx,1} = [];
        savedata.FtildeuTrans_Gbox{timeStepIdx,1} = [];
        savedata.eFtilde{timeStepIdx,1} = [];
    end
end

% particular solution
if isU
    savedata.G_PU_zero{timeStepIdx,1} = set.G_PU_zero;
    savedata.G_PU_infty{timeStepIdx,1} = set.G_PU_infty;
    savedata.PU_A_sum_error{timeStepIdx,1} = set.PU_A_sum_error;
else
    savedata.G_PU_zero{timeStepIdx,1} = [];
    savedata.G_PU_infty{timeStepIdx,1} = [];
    savedata.PU_A_sum_error{timeStepIdx,1} = [];
end

end

function [cnt,set,e,ebar,expmat,timeStep] = aux_ReadWriteLoop(k,cnt,...
    fullcomp,isU,set,e,ebar,expmat,timeStep,idx,readWrite)
% write to / read from cnti struct

if strcmp(readWrite,'read')
   
    % overwrite values for correct propagation
    timeStep = cnt.timeStep(idx);
    expmat.Deltatk = cnt.Deltatk{idx};
    
    e.acc = cnt.acc(idx);
    e.nonacc = cnt.nonacc(idx);
    ebar.acc(k) = cnt.ebaracc(idx);
    ebar.red(k) = cnt.ebarred(idx);
    ebar.rem(k) = cnt.ebarrem(idx);
    
    if fullcomp
        set.boxFstartset_center = cnt.boxFstartset_center{idx};
        set.boxFstartset_Gbox = cnt.boxFstartset_Gbox{idx};
        set.FtildeuTrans_center = cnt.FtildeuTrans_center{idx};
        set.FtildeuTrans_Gbox = cnt.FtildeuTrans_Gbox{idx};
    end
    if isU
        set.G_PU_zero = cnt.G_PU_zero{idx};
        set.G_PU_infty = cnt.G_PU_infty{idx};
    end
    
elseif strcmp(readWrite,'write')
    
    % save values for next cnt
    cnt.timeStep(idx,1) = timeStep;
    cnt.Deltatk{idx,1} = expmat.Deltatk;
    
    cnt.acc(idx,1) = e.acc;
    cnt.nonacc(idx,1) = e.nonacc;
    cnt.ebaracc(idx,1) = ebar.acc(k);
    cnt.ebarred(idx,1) = ebar.red(k);
    cnt.ebarrem(idx,1) = ebar.rem(k);
    
    if fullcomp
        cnt.boxFstartset_center{idx,1} = set.boxFstartset_center;
        cnt.boxFstartset_Gbox{idx,1} = set.boxFstartset_Gbox;
        cnt.FtildeuTrans_center{idx,1} = set.FtildeuTrans_center;
        cnt.FtildeuTrans_Gbox{idx,1} = set.FtildeuTrans_Gbox;
    end
    if isU
        cnt.G_PU_zero{idx,1} = set.G_PU_zero;
        cnt.G_PU_infty{idx,1} = set.G_PU_infty;
    end
    
end

end

% plot/debug functions to facilitate bug fixing
function debugdata = aux_debug_debugdata(k,cnt,fullcomp,debugdata,e,ebar,coeff,timeStep)
% save all current values for inspection in case bugs occur

debugdata.timeStep{k,1}(cnt,1) = timeStep;

% if fullcomp
% debugdata.coeff.linComb{k,1}.b(cnt,1) = coeff.linComb.b;
% debugdata.coeff.PU_tauk{k,1}.b(cnt,1) = coeff.PU_tauk.b;
% debugdata.coeff.F{k,1}.a(cnt,1) = coeff.F.a;
% debugdata.coeff.F{k,1}.b(cnt,1) = coeff.F.b;
% debugdata.coeff.Ftilde{k,1}.a(cnt,1) = coeff.Ftilde.a;
% debugdata.coeff.Ftilde{k,1}.b(cnt,1) = coeff.Ftilde.b;
% end
% debugdata.coeff.PU_tkplus1{k,1}.a(cnt,1) = coeff.PU_tkplus1.a;
% debugdata.coeff.PU_tkplus1{k,1}.b(cnt,1) = coeff.PU_tkplus1.b;

debugdata.e.nonacc{k,1}(cnt,1) = e.nonacc;
debugdata.e.linComb{k,1}(cnt,1) = e.linComb;
debugdata.e.PU_tauk{k,1}(cnt,1) = e.PU_tauk;
debugdata.e.F{k,1}(cnt,1) = e.F;
debugdata.e.Ftilde{k,1}(cnt,1) = e.Ftilde;
debugdata.e.acc{k,1}(cnt,1) = e.acc;
debugdata.e.red{k,1}(cnt,1) = e.red;
debugdata.e.PU_tkplus1{k,1}(cnt,1) = e.PU_tkplus1;

debugdata.ebar.acc{k,1}(cnt,1) = ebar.acc(k);
debugdata.ebar.red{k,1}(cnt,1) = ebar.red(k);
debugdata.ebar.rem{k,1}(cnt,1) = ebar.rem(k);
debugdata.ebar.remacc{k,1}(cnt,1) = ebar.remacc(k);
        
end

function debugdata = aux_debug_debugdata_sets(k,fullcomp,debugdata,set,expmat)
% save sets and matrices used for set propagation in each step

if fullcomp(k)
debugdata.set.startset{k,1} = set.startset;
debugdata.set.Hstartp{k,1} = set.Hstartp;
debugdata.set.Hstarti{k,1} = set.Hstarti;
end
debugdata.set.Pu{k,1} = set.Pu;
debugdata.set.Putotal{k,1} = set.Putotal;
debugdata.set.PU{k,1} = set.PU;
debugdata.set.PUtotal{k,1} = set.PUtotal;

debugdata.expmat.tk{k,1} = expmat.tk;
debugdata.expmat.tkplus1{k,1} = expmat.tkplus1;
debugdata.expmat.Deltatk{k,1} = expmat.Deltatk;

end

function aux_plot_Rerror(isdebug,fullcomp,tVec,Rout_error,Rout_tp_error,e,tu)
% plot curves of over-approximation error contained in time-point and
% time-interval reachable sets, along with maximum admissible error and
% switches of input vector (to possibly explain some jumps)

% isdebug = true;
if ~isdebug
    return;
end

% generate step function (repeat the inner values)
cumsumtVec_tp = cumsum(tVec);
cumsumtVec = [0;repelem(cumsumtVec_tp(1:end-1),2);cumsumtVec_tp(end)];

Rout_error_step = repelem(Rout_error(fullcomp),2);

figure; hold on; box on;
axis([0,cumsumtVec_tp(end),0,1.2*e.emax]);
% plot errors in time-interval and time-point reachable sets
h_ti = plot(cumsumtVec,Rout_error_step,'Color',colorblind('b'),'LineWidth',1.5);
h_tp = scatter(cumsumtVec,Rout_tp_error(fullcomp),8,colorblind('b'),'filled');
% maximum error
hBound = plot([0;cumsumtVec_tp(end)],[e.emax;e.emax],'Color',colorblind('r'),'LineWidth',2);
% plot vertical lines for changes in input vector (if given)
if length(tu) > 1
    for u=2:length(tu)
        h_u = plot([tu(u);tu(u)],[0;1.2*e.emax],'Color',colorblind('gray'));
    end
    legend([hBound,h_ti,h_tp,h_u],'$\varepsilon_{\max}$','$\varepsilon(\mathcal{R}(\tau_k))$',...
        '$\varepsilon(\mathcal{R}(t_k))$','$u_k$','interpreter','latex');
else
    legend([hBound,h_ti,h_tp],'$\varepsilon_{\max}$','$\varepsilon(\mathcal{R}(\tau_k))$',...
        '$\varepsilon(\mathcal{R}(t_k))$','interpreter','latex');
end
% labels
xlabel('$t$','interpreter','latex');
ylabel('$\varepsilon$','interpreter','latex');
% matlab2tikz('fig_errorbounds.tikz');
close;

end

function aux_plot_errorcurve(isdebug,fullcomp,e,ebar,tVec,debugdata)
% plots the committed error and the corresponding bounds over time

% isdebug = true;
if ~isdebug
    return;
end

steps = length(fullcomp);
% rewrite into readable arrays
ebarnonacc = ebar.rem;
enonacc = zeros(steps,1); eacc = zeros(steps,1); ered = zeros(steps,1);
for i=1:steps
    enonacc(i,1) = debugdata.e.nonacc{i}(end);
    eacc(i,1) = debugdata.e.acc{i}(end);
    ered(i,1) = debugdata.e.red{i}(end);
end

% cumulative sum to yield bound for accumulating error
% note: subtract error committed in the step!
ebaracctotal = e.acctotal + ebar.acc - eacc;
ebarredtotal = e.redtotal + ebar.red - ered;

% compute total committed error (time-interval solution)
totalerror = enonacc + [0;e.acctotal(1:end-1)] + e.redtotal;
totalerror = totalerror(fullcomp);

% generate time vector
cumsumtVec_tp = cumsum(tVec);
cumsumtVec = [0;repelem(cumsumtVec_tp(1:end-1),2);cumsumtVec_tp(end)];

% start and end of time horizon
tStart = 0;
tFinal = cumsumtVec_tp(end);

rows = 2; cols = 3; lw = 0.5;
figure;

% 1. total error
subplot(rows,cols,1); hold on; box on;
title("total error (ti)");
axis([tStart,tFinal,0,1.1*e.emax]);
% committed error
plot(cumsumtVec,repelem(totalerror,2),'Color',colorblind('b'),'LineWidth',lw);
% maximum error (=bound)
plot([0;tFinal],[e.emax;e.emax],'Color',colorblind('r'),'LineWidth',lw);


% 2. accumulating error (single step)
subplot(rows,cols,2); hold on; box on;
title("eacc, single step");
% committed error
plot(cumsumtVec,repelem(eacc,2),'Color',colorblind('b'),'LineWidth',lw);
% error bound
plot(cumsumtVec,repelem(ebar.acc,2),'Color',colorblind('y'),'LineWidth',lw);


% 3. accumulating error (accumulated)
subplot(rows,cols,3); hold on; box on;
title("eacc, accumulated");
axis([tStart,tFinal,0,1.1*e.emax]);
% committed error
plot(cumsumtVec,repelem(e.acctotal,2),'Color',colorblind('b'),'LineWidth',lw);
% error bound
plot(cumsumtVec,repelem(ebaracctotal,2),'Color',colorblind('y'),'LineWidth',lw);
% maximum error
plot([0;tFinal],[e.emax;e.emax],'Color',colorblind('r'),'LineWidth',lw);


% 4. non-accumulating error
subplot(rows,cols,4); hold on; box on;
title("enonacc");
axis([tStart,tFinal,0,1.1*e.emax]);
% committed error
plot(cumsumtVec,repelem(enonacc(fullcomp),2),'Color',colorblind('b'),'LineWidth',lw);
% error bound
plot(cumsumtVec,repelem(ebarnonacc(fullcomp),2),'Color',colorblind('y'),'LineWidth',lw);
% maximum error
plot([0;tFinal],[e.emax;e.emax],'Color',colorblind('r'),'LineWidth',lw);


% 5. reduction error (single step)
subplot(rows,cols,5); hold on; box on;
title("ered, single step");
% committed error
plot(cumsumtVec,repelem(ered,2),'Color',colorblind('b'),'LineWidth',lw);
% error bound
plot(cumsumtVec,repelem(ebar.red,2),'Color',colorblind('y'),'LineWidth',lw);
% error bound curve
plot(ebar.red_t,ebar.red_e,'Color',colorblind('r'),'LineWidth',lw);


% 6. reduction error (accumulated)
subplot(rows,cols,6); hold on; box on;
title("ered, accumulated");
axis([tStart,tFinal,0,1.1*e.emax]);
% committed error
plot(cumsumtVec,repelem(e.redtotal,2),'Color',colorblind('b'),'LineWidth',lw);
% error bound
plot(cumsumtVec,repelem(ebarredtotal,2),'Color',colorblind('y'),'LineWidth',lw);
% error bound curve
plot(ebar.red_t,ebar.red_e,'Color',colorblind('r'),'LineWidth',lw);
% maximum error
plot([0;tFinal],[e.emax;e.emax],'Color',colorblind('r'),'LineWidth',lw);

close;

end

function aux_plot_e2ebar(isdebug,fullcomp,debugdata,timeStepIdxs)
% plot to what percentage the bounds are fulfilled in each step

% isdebug = true;
if ~isdebug
    return;
end

% pre-process plotting data
steps = length(fullcomp);
enonacc_sat = zeros(steps,1);
eacc_sat = zeros(steps,1);
% retrieve data from saved data
for k=1:steps
    if fullcomp(k)
        enonacc_sat(k) = debugdata.e.nonacc{k}(timeStepIdxs(k)) / debugdata.ebar.rem{k}(timeStepIdxs(k));
    end
    eacc_sat(k) = debugdata.e.acc{k}(timeStepIdxs(k)) / debugdata.ebar.acc{k}(timeStepIdxs(k));
end
% get percentages
enonacc_sat = 100 * enonacc_sat;
eacc_sat = 100 * eacc_sat;

figure; hold on; grid on;
h_enonacc = plot(enonacc_sat,'Color',colorblind('b'));
h_eacc = plot(eacc_sat,'Color',colorblind('r'));
axis([0,steps+1,-1,101]);
xlabel('Step');
ylabel('Fulfillment percentage');
legend([h_enonacc,h_eacc],'enonacc','eacc',...
    'Location','southeast');
close;

end

function aux_print_e2ebar(isdebug,k,cnt,isU,fullcomp,e,ebar,timeStep)

% isdebug = true;
if ~isdebug
    return;
end

% print time step size
disp("(" + k + "/" + cnt + ") - Delta t: " + timeStep);
% print fulfillment of error bound in percent
if fullcomp
disp("(" + k + "/" + cnt + ") - erem: " + ...
    sprintf('%0.2f',round(100 * (e.nonacc + e.acc) / ebar.rem(k),3)) + "%");
end
if isU
disp("(" + k + "/" + cnt + ") - eacc: " + ...
    sprintf('%0.2f',round(100 * e.acc / ebar.acc(k),3)) + "%");
end

end

% ------------------------------ END OF CODE ------------------------------
