function [Rti,Rtp,options] = linReach_adaptive(obj,options,Rstart)
% linReach_adaptive - computes the reachable set after linearization;
%    based on automated tuning from [3] applied to algorithms [1,2]
%
% Syntax:
%    [Rti,Rtp,options] = linReach_adaptive(obj,options,Rstart)
%
% Inputs:
%    obj     - nonlinearSys or nonlinearParamSys object
%    options - struct with algorithm settings
%    Rstart  - reachable set (time point of current time)
%
% Outputs:
%    Rti - reachable set for time interval
%    Rtp - reachable set for time point
%    options - struct with algorithm settings
%
% References:
%   [1] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%       uncertain parameters using conservative linearization"
%   [2] M. Althoff et al. "Reachability analysis of nonlinear systems using
%       conservative polynomialization and non-convex sets"
%   [3] M. Wetzlinger et al. "Automated parameter tuning for reachability
%        analysis of nonlinear systems", HSCC 2021.
%   [4] M. Wetzlinger et al. "Adaptive reachability algorithms for
%        nonlinear systems using abstraction error analysis", NAHS 2022.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       11-December-2020
% Last update:   15-January-2021
%                30-June-2021
%                10-November-2023 (MW, improved estimate for finitehorizon)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set linearization error
error_adm = options.error_adm_horizon;

% time step size + set abstraction order ----------------------------------
% finite horizon depending on varphi, init at start, control below
lastStep = false; veryfirststep = false; timeStepequalHorizon = false;

% first step
if options.i == 1
    veryfirststep = true;
    % initial guess for time step size / finite horizon
    options.timeStep = (options.tFinal - options.tStart) * 0.01;
    finitehorizon = options.timeStep;

    options.tt_err = [];

    % init reduction indices (only used for options.tensorOrder = 3)
    if strcmp(options.alg,'lin')
        options.gredIdx.Rhomti = {};
        options.gredIdx.Rhomtp = {};
        options.gredIdx.Rpar = {};
        options.gredIdx.Rred = {};
        options.gredIdx.VerrorDyn = {};
    end

% last step
elseif options.timeStep > options.tFinal - options.t
    lastStep = true;
    options.timeStep = options.tFinal - options.t;
    finitehorizon = options.timeStep;
    
% non-start/end step
elseif options.i > 1

    % take last one and estimate so that varphi ~ zetaphi
    % new estimation (better than line 1 from [4, Alg. 2])
    finitehorizon = options.finitehorizon(options.i-1) ...
        * (1 + options.varphi(options.i-1) - options.zetaphi(options.minorder+1));

    % finitehorizon is capped by remaining time
    min([options.tFinal - options.t, finitehorizon]);
    
    % addition for ARCH'23: artificially enlarge time step if current
    % finitehorizon almost same as last time step size
    if withinTol(finitehorizon,options.stepsize(options.i-1),1e-3)
        finitehorizon = options.finitehorizon(options.i-1) * 1.1;
    end

    % init time step size for current step as finitehorizon
    options.timeStep = finitehorizon;
    
    assert(options.timeStep > 0,'Tuning error.. report to devs');
    
    if options.i == 2 && strcmp(options.alg,'poly') && options.polyGain
        % gain of one order in poly-algorithm is lost after one step
        % (exact reasons unclear as of now)
        options.orders = options.orders - 1;
        options.minorder = min(options.orders);
    end
end

% check if Rstart is degenerate (each dimension): order may change
zeroWidthDim = withinTol(sum(abs(generators(zonotope(Rstart))),2),0);

% iteration counter run:
% - 0: only for first step in lin, otherwise:
% - 1: compute finitehorizon (use last deltat for varphi),
% - 2: compute optimal deltat
options.run = 0;

% first step: adaptive tensor order
if veryfirststep && strcmp(options.alg,'lin')
    % compute L for start set to decide abstraction order
    options = aux_initStepTensorOrder(obj,options,Rstart);
end
% thistO = options.tensorOrder;
% nexttO = options.tensorOrder; % init
options.run = options.run + 1;
% -------------------------------------------------------------------------

% log
abscount = 0;

% loop for time step adaptation
while true
    
    % linearize the nonlinear system
    [obj,linSys,linOptions] = linearize(obj,options,Rstart);
    
    if options.run == 1
        % scaling over delta t, required for optimization function
        zetaP = exp(trace(linSys.A*options.timeStep));
    end
    
    % translate Rstartset by linearization point
    Rdelta = Rstart + (-obj.linError.p.x);

    % compute reachable set of the linearized system
    % using @linSys > initReach
    try
        [Rlin,options,linOptions] = initReach_adaptive(linSys,Rdelta,options,linOptions);
    catch ME % if time step size is too big, then Inf/NaN in expmtie
        if strcmp(ME.identifier,'expmtie:notconverging')
            options.timeStep = options.timeStep * 0.5;
            finitehorizon = options.timeStep;
            continue
        else
            rethrow(ME);
        end
    end

    % compute static error out of loop below
    if strcmp(options.alg,'poly') && options.run == 1
        [H,Zdelta,errorStat,T,ind3,Zdelta3] = ...
        	precompStatError_adaptive(obj,options,Rdelta);
    end
    
    % compute reachable set of the abstracted system including the
    % abstraction error using the selected algorithm

    % loop until the actual abstraction error is smaller than
    % the estimated abstraction error
    Rlintp = Rlin.tp; Rlinti = Rlin.ti;
    
    perfIndCounter = 1; perfInds = []; options.Lconverged = false;
    while true
        % increment counter how often abstraction error computed
        abscount = abscount + 1;
        
        % set resulting from (estimated!) abstraction error
        if any(error_adm)
            errG = diag(error_adm);
            Verror = zonotope([0*error_adm,errG(:,any(errG,1))]);
            % compute abstraction error as input solution
            [RallError,options] = errorSolution_adaptive(linSys,options,Verror);
        else
            % only used in options.run == 1, but important for adaptive
            RallError = zonotope(zeros(obj.dim,1));
        end
        
        try
            % compute overall reachable set including abstraction error
            if strcmp(options.alg,'lin')
                Rmax = Rlinti + RallError;
                Rdiff = zonotope(zeros(length(center(RallError)),1));
                Rtemp = Rdiff + RallError;

                % compute abstraction error
                % (dependency on settings resolved inside)
                [VerrorDyn, VerrorStat, trueError, options] = ...
                    abstractionError_adaptive(obj,options,Rmax,Rtemp);
                
                
            elseif strcmp(options.alg,'poly')
                Rdiff = deltaReach_adaptive(linSys,Rdelta,linOptions); % no red
                try   Rmax = Rdelta+zonotope(Rdiff)+RallError;
                catch Rmax = Rdelta+Rdiff+RallError; 
                end
                % compute abstraction error
                % (dependency on settings resolved inside)
                [VerrorDyn, VerrorStat, trueError, options] = ...
                    abstractionError_adaptive(obj,options,Rmax,Rdiff+RallError,...
                    H,Zdelta,errorStat,T,ind3,Zdelta3);

            else
                error("Unsupported options.alg");
            end
        catch ME
            if strcmp(ME.identifier,'reach:setoutofdomain')
                options.Lconverged = false; break
            else
                rethrow(ME);
            end
        end
        
        % compare linearization error with the maximum admissible error
        perfIndCurr = max(trueError ./ error_adm);
        if perfIndCurr <= 1 || ~any(trueError)
            % just for information
            perfInds(perfIndCounter) = perfIndCurr;
            options.Lconverged = true;
            break
        elseif perfIndCounter > 1
            perfInds(perfIndCounter) = perfIndCurr;
            if perfIndCounter > 2 && perfInds(perfIndCounter) > perfInds(perfIndCounter-1)
                options.Lconverged = false; break
            end
        end
        
        % increase admissible abstraction error for next iteration
        error_adm = 1.1 * trueError;
        perfIndCounter = perfIndCounter + 1;
    end
    
    % if either L was out of domain or did not converge
    if ~options.Lconverged
        % half time step, reset linearization error
        options.timeStep = options.timeStep * 0.5;
        finitehorizon = options.timeStep;
        error_adm = zeros(obj.dim,1); % options.error_adm_horizon
        continue
    end
    % ... now containment of L ensured
    
    % if rank of Rstart changes... re-compute orders
    if veryfirststep || ~all(zeroWidthDim == options.zeroWidthDim)
        options.zeroWidthDim = zeroWidthDim;
        options = aux_getPowers(obj,options,linSys,zeroWidthDim,VerrorDyn);
        % old version
%         options = aux_getPowers_old(linSys,options,VerrorDyn);
    end
    
    % compute set of abstraction errors
    Rerror = errorSolution_adaptive(linSys,options,VerrorDyn,VerrorStat);
    
    % measure abstraction error
    if isa(Rerror,'polyZonotope')
        abstrerr = sum(abs(Rerror.GI),2)';
        % Here, we take the independent generators, as they represent
        % relatively truthfully the part of the abstraction error
        % corresponding to the V^{\Delta} (see [2] Section 4.2). However,
        % there is also a part of the static error remaining in this part,
        % which contaminates the result if Rerror is a polynomial zonotope. 
        % This leads to an order decrease of 1, which happens to be 
        % beneficial for our approach, as it ensures the stability of the 
        % time step choosing (i.e., the chosen time-step is not too small).
    else
        abstrerr = sum(abs(generators(Rerror)),2)';
    end
    
    % two runs in each step
    if options.run == 1 && ~lastStep
        % ... run using finite horizon
        
        % first step more cumbersome as there is no previous knowledge
        if options.i == 1
            
            if veryfirststep
                veryfirststep = false;
                % save values from finite horizon
                abstrerr_h = abstrerr;
                Rerror_h = Rerror;
                Rti_h = Rlinti; Rtp_h = Rlintp;
                linx_h = obj.linError.p.x;
                zetaP_h = zetaP;
                % decrease time step size and linearization error
                options.timeStep = options.decrFactor * options.timeStep;
                error_adm = options.decrFactor * trueError;
                continue;
            end
            
            % take worst-case current gain (use same dimensions as for lin)
            % to avoid numerical issues (i.e., lower orders jumping to
            % "impossible" values for varphi_i due to other effects), take
            % only dimensions with lowest order
            temp = abstrerr(options.orders == options.minorder) ./ ...
                abstrerr_h(options.orders == options.minorder);
            varphi = max( temp(~isnan(temp)) );
            
            % check condition for varphi
            if varphi < options.zetaphi(options.minorder+1)
                % decrease time step size
                finitehorizon = options.timeStep;
                options.timeStep = options.decrFactor * options.timeStep;
                % reset estimate for linearization error
                error_adm = options.decrFactor * trueError;
                % save values from finite horizon
                abstrerr_h = abstrerr;
                Rerror_h = Rerror;
                Rti_h = Rlinti; Rtp_h = Rlintp;
                linx_h = obj.linError.p.x;
                zetaP_h = zetaP;
                continue;
            end
            
            options.varphi(options.i,1) = varphi;
            options.finitehorizon(options.i,1) = finitehorizon;
            % estimate near-optimal delta t by optimization function
            [options.timeStep,~] = aux_optimaldeltat(Rstart,Rerror_h,...
                finitehorizon,varphi,zetaP_h,options);
            % save abstraction error, reset to zero for optimal time step size
            options.error_adm_horizon = trueError;
            error_adm = zeros(obj.dim,1);
            
        % non-start/end step
        elseif options.i > 1
            % save values from finite horizon
            abstrerr_h = abstrerr;
            Rerror_h = Rerror;
            Rti_h = Rlinti; Rtp_h = Rlintp;
            linx_h = obj.linError.p.x;
            % solve optimization problem for time step size
            [options.timeStep,~] = aux_optimaldeltat(Rstart,Rerror_h,...
                finitehorizon,options.varphi(options.i-1),zetaP,options);
            % set time step to remaining time horizon if too large
            options.timeStep = min([options.timeStep, options.tFinal - options.t]);
            % save linearization error for next step
            options.error_adm_horizon = trueError;
            % set linearization error to values for previous optimal delta t
            error_adm = options.error_adm_Deltatopt;
        end
        
    elseif options.run == 2 || lastStep
        % ... run using tuned time step size (or scaled time step size
        % for varphi if timeStepEqualHorizon is true, or last step)
        if timeStepequalHorizon
            temp = abstrerr(options.orders == options.minorder) ./ ...
                abstrerr_h(options.orders == options.minorder);
            options.varphi(options.i,1) = max( temp(~isnan(temp)) );
        elseif ~lastStep
            options.varphi(options.i,1) = aux_varphiest(finitehorizon,...
                options.timeStep,Rerror_h,Rerror,options.decrFactor,...
                options.orders,options.minorder);
        end
        % save finite horizon and linearization error for next step
        options.finitehorizon(options.i,1) = finitehorizon;
        options.error_adm_Deltatopt = trueError;
        
        % predict tensor order of next step
        if strcmp(options.alg,'lin')
            if options.i == 1
                % save time step and abstrerr
                options.kappa_deltat = options.timeStep;
                options.kappa_abstrerr = vecnorm(abstrerr);
            elseif ~lastStep
                % compensate for the difference between the current time
                % step size and stored time step size by varphi
                options = aux_nextStepTensorOrder(obj,options,...
                    vecnorm(abstrerr),linSys,Rmax,Rtemp);
            end
        end
        
        % second run finished -> exit
        break
    end
    
    % update run counter
    options.run = options.run + 1;
    
    if options.timeStep == finitehorizon
        % required Rerror and Rlin already there, but compute
        % timeStep = decrFactor*finitehorizon for usable varphi
        % (otherwise two points required for varphi are the same point)
        timeStepequalHorizon = true;
        options.timeStep = options.timeStep * options.decrFactor;
    end 
    
end

% use correct Rerror (if optimal time step size is equal to finitehorizon)
if timeStepequalHorizon
    % reset time step size for correct time increment in loop outside
    options.timeStep = finitehorizon;
    % choose sets computed by horizon (last run: finitehorizon*decrFactor)
    Rti = Rti_h + linx_h; Rtp = Rtp_h + linx_h;
    Rerror = Rerror_h;
else
    % translate reachable sets by linearization point
    Rti = Rlinti + obj.linError.p.x; Rtp = Rlintp + obj.linError.p.x;
end

% add the set of abstraction errors to the reachable sets
if isa(Rerror,'polyZonotope')
    Rti = exactPlus(Rti,Rerror);
    Rtp = exactPlus(Rtp,Rerror);
else
    Rti = Rti + Rerror; Rtp = Rtp + Rerror;
end

% time step equal to finite horizon
options.timeStepequalHorizon(options.i,1) = timeStepequalHorizon;
% save counter of how often Lagrange remainder was computed
options.abscount(options.i,1) = abscount;
% save time step size
options.stepsize(options.i,1) = options.timeStep;

end


% Auxiliary functions -----------------------------------------------------

% adaptation of tensorOrder: init step
function options = aux_initStepTensorOrder(obj,options,Rstartset)
% computes L for the start set of the first step to determine tensorOrder

% special handling for initial set as just a single point
if isempty(generators(Rstartset))
    options.tensorOrder = 2;
    return;
end

obj.linError.p.u = center(options.U);
obj.linError.p.x = center(Rstartset);
Rdelta = Rstartset + (-obj.linError.p.x);

% compute L for both kappas
options.tensorOrder = 2;
[~,~,L0_2,options] = abstractionError_adaptive(obj,options,Rdelta);
options.tensorOrder = 3;
[~,~,L0_3,options] = abstractionError_adaptive(obj,options,Rdelta);

% cut linear dimensions
L0_2lin = L0_2(L0_2 ~= 0);
L0_3lin = L0_3(L0_3 ~= 0);

% compare Lagrange remainder to decide which tensorOrder to use
if all( L0_3lin ./ L0_2lin > options.zetaK )
    options.tensorOrder = 2;
else
    options.tensorOrder = 3;
end


end

% adaptation of tensorOrder: every step for next step
function options = aux_nextStepTensorOrder(obj,options,abstrerr,linSys,Rmax,Rtemp)

% time step sizes
timeStep = options.timeStep;
timeStep_prev = options.kappa_deltat;

% compute varphi for given time step sizes
varphi_lim = options.decrFactor^(options.minorder+1);
varphitotal = aux_varphitotal(options.varphi(options.i),varphi_lim,...
    timeStep_prev,timeStep,options.decrFactor);

% abstraction error
abstrerr_prev = options.kappa_abstrerr;
% compute estimate of abstraction error
abstrerr_est = varphitotal * abstrerr;

% compare to last saved abstraction error (computed using same kappa)
if options.tensorOrder == 2 
    % original criterion:
    % only check other tensorOrder if abstraction error has grown
%     if abstrerr_est / abstrerr_prev - 1 > 1 - options.zetaK
    % adapted criterion:
    if abs(abstrerr_est / abstrerr_prev - 1) > 1 - options.zetaK
        abstrerr_2 = abstrerr;
        % compute abstraction error for other tensor order
        options.tensorOrder = 3;
    else
        return;
    end
    
    % compute set of abstraction errors
    [VerrorDyn,~,err] = abstractionError_adaptive(obj,options,Rmax,Rtemp);
    Rerror = errorSolution_adaptive(linSys,options,VerrorDyn);
    % simplify to scalar value
    abstrerr_3 = vecnorm(sum(abs(generators(Rerror)),2));
    
    % compare results
    if ~all( abstrerr_3 ./ abstrerr_2 > options.zetaK )
        options.kappa_deltat = options.timeStep;
        options.kappa_abstrerr = abstrerr_3;
        options.error_adm_Deltatopt = err;
        options.error_adm_horizon = zeros(obj.dim,1);
    else
        options.tensorOrder = 2;
    end
    
elseif options.tensorOrder == 3
    % original criterion:
    % only check other tensorOrder if abstraction error has shrunk
%     if 1 - abstrerr_est / abstrerr_prev < 1 - options.zetaK
    % adapted criterion:
    if abs(abstrerr_est / abstrerr_prev - 1) > 1 - options.zetaK
        abstrerr_3 = abstrerr;
        % compute other tensor order and check again
        options.tensorOrder = 2;
    else
        return;
    end
    
    % compute set of abstraction errors
    [VerrorDyn,~,err] = abstractionError_adaptive(obj,options,Rmax,Rtemp);
    Rerror = errorSolution_adaptive(linSys,options,VerrorDyn);
    % simplify to scalar value
    abstrerr_2 = vecnorm(sum(abs(generators(Rerror)),2));
    
    % compare results
    if all( abstrerr_3 ./ abstrerr_2 > options.zetaK )
        options.kappa_deltat = options.timeStep;
        options.kappa_abstrerr = abstrerr_2;
        options.error_adm_Deltatopt = err;
        options.error_adm_horizon = zeros(obj.dim,1);
    else
        options.tensorOrder = 3;
    end
end

end

% compute total scaling for two time step sizes given some starting varphi
function varphitotal = aux_varphitotal(varphi,varphi_lim,timeStep_prev,timeStep_curr,decrFactor)

% number of times current time step size needs to be scaled in order to
% reach previous time step size
nrScalings = log(timeStep_prev / timeStep_curr) / log(decrFactor);

% since varphi is only valid for fixed decrFactor / at current time step,
% we have to compute the remaining list of varphis by linear interpolation
varphis = varphi;

% current time step size smaller than previous one
if nrScalings > 0
    
    nrScalings_low = floor(nrScalings + 10*eps);
    nrScalings_high = nrScalings_low + 1;
    
    % linearly interpolate between (0,varphi_lim) and (timeStep,varphi)
    for i=1:nrScalings_low
        varphis(i+1,1) = varphi + (varphi_lim - varphi) * ...
            (timeStep_curr - decrFactor^i * timeStep_curr) / (timeStep_curr - 0);
    end
    
    % compute remaining partial decrFactor
    varphis(end) = varphis(end) + (1 - varphis(end)) * ...
        (timeStep_prev - timeStep_curr * decrFactor^nrScalings_high) / ...
        (timeStep_curr * (decrFactor^nrScalings_low - decrFactor^nrScalings_high));
    
% current time step size larger than previous one
else
    
    nrScalings_low = ceil(nrScalings - 10*eps);
    nrScalings_high = nrScalings_low - 1;
    
    % linearly extrapolate the line (0,varphi_lim) -- (timeStep,varphi)
    for i=-1:-1:nrScalings_low
        varphis(end+1,1) = varphi + (varphi_lim - varphi) * ...
            (timeStep_curr - decrFactor^i * timeStep_curr) / (timeStep_curr - 0);
    end
    
    % compute remaining partial decrFactor
    varphis(end) = varphis(end) + (1 - varphis(end)) * ...
    	(timeStep_prev - timeStep_curr * decrFactor^nrScalings_high) / ...
        (timeStep_curr * (decrFactor^nrScalings_low - decrFactor^nrScalings_high));
    
    varphis = 1 ./ varphis;
    
end

% combine varphis
varphitotal = prod(varphis);

end

% optimization function for optimal time step size
function [deltatest,kprimeest] = aux_optimaldeltat(Rt,Rerr,deltat,varphimin,zetaP,opt)

% read from options struct
mu = opt.decrFactor;
if strcmp(opt.alg,'lin')
    dHused = 0.5;
elseif strcmp(opt.alg,'poly')
    % since poly over-approximation is a lot more conservative
    dHused = 0.3;
end
zetaZ = opt.redFactor * dHused; % better representation of actual reduction

% kprimemax = log(deltatmin/deltat) / log(mu); % if deltatmin computed
kprimemax = ceil(-log(100) / log(mu)); % max 100 steps instead of 1 (10000 total steps max)
kprime = 0:round(kprimemax); % just to catch rounding errors
k = mu .^ (-kprime);
deltats = deltat * mu.^kprime; % consecutively scaled down by mu
% number of "full" steps, i.e., steps where deltats(i) fully used (all but last)
floork = floor(k);

if isa(Rt,'zonotope') && isa(Rerr,'zonotope')
    rR = vecnorm(sum(abs(generators(Rt)),2),2);
    rerr1 = vecnorm(sum(abs(generators(Rerr)),2),2);
else
    rR = vecnorm(rad(interval(Rt)));
    rerr1 = vecnorm(sum(abs(Rerr.GI),2),2);
end

% model varphi by linear interpolation over the gain
% ...bounds: 1 (at deltatmin) to measured gain (at deltat)
varphimax = mu; varphi_h = (varphimax - varphimin);
varphi = (varphimin + (deltats(1) - deltats)/deltats(1) * varphi_h) / mu;
varphiprod = cumprod(varphi);

% optimization by minimization of set size after finite horizon
sumallbutlast = zeros(1,length(floork));
for i=1:length(floork)
    % evaluate sum for each floork
    firstfactor = (1+2*zetaZ) .^ (k(i) + 1 - (1:floork(i)));
    secondfactor = zetaP .^ ( 1 - (1:floork(i)) ./ k(i));
    sumallbutlast(i) = sum(firstfactor .* secondfactor);
end

% zetaZ = opt.redFactor; % adapt for high-dimensional systems
objfuncset = rR * (1+2*zetaZ) .^ k * zetaP + rerr1 ./ k .* varphiprod .* ...
    ( sumallbutlast + (1+zetaZ).^(k - kprime) .* (k - floork) );
[~,bestIdxnew] = min(objfuncset);
deltatest = deltats(bestIdxnew);
kprimeest = bestIdxnew - 1;

% rabskest = rerr1 / k(bestIdxnew) * varphiprod(bestIdxnew);

end

% estimation of order for given finite horizon
function varphi = aux_varphiest(horizon,deltat,Rerr_h,Rerr_deltat,...
    decrFactor,orders,minorder)

% compute estimate for set sizes
if isa(Rerr_h,'zonotope') && isa(Rerr_deltat,'zonotope')
    % read out generator matrices
    G_Rerr_h = generators(Rerr_h);
    G_Rerr_deltat = generators(Rerr_deltat);
    % only least-order directions
    rerr1 = vecnorm(sum(abs(G_Rerr_h(orders == minorder,:)),2),2);
    rerrk = vecnorm(sum(abs(G_Rerr_deltat(orders == minorder,:)),2),2);
    % radius
%     rerr1 = vecnorm(sum(abs(generators(Rerr_h)),2),2);
%     rerrk = vecnorm(sum(abs(generators(Rerr_deltat)),2),2);
else
    rerr1 = vecnorm(sum(abs(Rerr_h.GI),2),2);
    rerrk = vecnorm(sum(abs(Rerr_deltat.GI),2),2);
end

% sanity check: rerr1 and rerrk are computed in the same time step,
% therefore rerr1 has to be larger than rerrk unless there are some
% artifacts or we are in a theoretically unexplored region of Rerror
assert(rerr1 > rerrk,'Check abstraction errors');

% right-hand side ... total scaling by varphi
rhs = rerrk / rerr1;

% limit varphi
varphi_lim = decrFactor^(minorder+1); % if power always 1, former results

% use bisection to find varphi
% upper bound: varphi = limit value at deltat -> 0
varphi_up = decrFactor^(minorder+1);
% compute accumulated scaling
varphitotal = aux_varphitotal(varphi_up,varphi_lim,deltat,horizon,decrFactor);
% compute residual of prod of varphi (assuming varphi) and rhs
% residual_up = varphitotal - rhs;

% lower bound: varphi = 0
varphi_low = 0;
% compute accumulated scaling
varphitotal = aux_varphitotal(varphi_low,varphi_lim,deltat,horizon,decrFactor);
% compute residual of prod of varphi (assuming varphi) and rhs
% residual_low = varphitotal - rhs;

% counter for number of iterations
cnt = 0;

while true
    
    % increment counter
    cnt = cnt + 1;
    
    % update varphi by bisection
    varphi(cnt,1) = varphi_low + 0.5 * (varphi_up - varphi_low);
    % old version
%     varphi(cnt,1) = varphi_low + (varphi_up - varphi_low) * ...
%         -residual_low / (residual_up - residual_low);

    % compute accumulated scaling
    varphitotal = aux_varphitotal(varphi(cnt),varphi_lim,deltat,horizon,decrFactor);

    % compute residual of prod of varphi (assuming varphi) and rhs
    residual(cnt,1) = varphitotal - rhs;
    if residual(cnt) < 0
%         residual_low = residual(cnt);
        varphi_low = varphi(cnt);
    else
%         residual_up = residual(cnt);
        varphi_up = varphi(cnt);
    end
    
    % break condition
    if cnt == 10000
        error("Bug in varphi estimation... report to devs");
    end
    if abs(residual(cnt)) < 1e-9 || ...
            (cnt > 1 && abs(varphi(cnt) - varphi(cnt-1)) < 1e-6)
        varphi = varphi(cnt);
        break;
    end

end

end

% order of abstraction error
function options = aux_getPowers(obj,options,linSys,zeroWidthDim,Vdyn)
% only called once after the very first converged L!

% hard-coded shortcut
if options.i == 1 && isfield(options,'orders')
    options.minorder = min(options.orders);
    return;
end

% propagation matrix
A = linSys.A;
n = obj.dim;
% m = obj.nrOfInputs;
    
% start set is just a point
if all(zeroWidthDim)
    sigma = 2*ones(n,1);

% start set is full-dimensional (non-degenerate)
elseif ~any(zeroWidthDim)
    sigma = zeros(n,1);

% start set is degenerate
elseif any(zeroWidthDim)

    % examine every row in G to obtain dependency on time step size;
    % this is done by going over all products z' * H{i} * z
    % and checking whether the resulting term converges to 0 and if so,
    % how the zeroWidthDim-dimensions linear or quadratic factors
    sigma = zeros(n,1); % init
    sigmafound = false(n,1);
    
    % evaluate Hessian...
    % assumptions:  1. no order over time
    %               2. no accidental zeros
    % ... hence, evaluation on linearization point sufficient
    % TODO: make this more robust if needed
    if options.isHessianConst
        % easier case: Hessian is constant
        Hess = options.hessianConst;
    else
        if strcmp(options.alg,'poly')
            Hess = obj.hessian(obj.linError.p.x,obj.linError.p.u);
        else
            % Hessian requires evaluation on intervals for x and u...
            % try eval on linearization point, set handle to correct file
            obj = setHessian(obj,'standard');
            Hess = obj.hessian(obj.linError.p.x,obj.linError.p.u);
        end
    end
    
    for i=1:n

        % take Hessian corresponding to current dimension
        H = Hess{i};
        
        % go over non-zeroWidthDims
        for j=1:n
            if ~zeroWidthDim(j) && H(j,j) ~= 0
                sigma(i) = 0; sigmafound(i) = true; break;
            end
        end
        % go to next dimension
        if sigmafound(i)
            continue;
        end

        % go over sum of off-diag entries
        sigma(i) = Inf;
        for j=1:n
            for jj=j+1:n
                if H(j,jj) + H(jj,j) ~= 0
                    sigma(i) = min([sigma(i), nnz(zeroWidthDim([j,jj]))]);
                    sigmafound(i) = true;
                end
            end
        end

        % sigma still not found (should not be the case)
        if ~sigmafound(i)
            sigma(i) = 2;
        end

    end
end

% poly: second-order dynamical error can go to 0 by O(t^1)
% ... only if L(3) is all-zero, otherwise dependencies...
% TODO: extend to arbitrary third-order tensors (no cases in current
% systems so can be done later)
if strcmp(options.alg,'poly')
    options.polyGain = false;
    if options.thirdOrderTensorempty
        sigma = sigma + 1;
        options.polyGain = true;
    end
end

% express G as function over sigma (implicitly)
Gzero = vecnorm(generators(Vdyn),1,2) == 0;

% compute q_i
% ... define G depending on sigmas (dependency on t!)
qi = Inf(n,1);
p = 0;
Apower = eye(n);
qip = [];
while true
    
    % loop over q_i
    for i=1:n
        % init q_i
        qip(i,p+1) = Inf;
        for j=1:n
            % loop over all multiplications
            if ~Gzero(j) && Apower(i,j) ~= 0
                qip(i,p+1) = min([qip(i,p+1),sigma(j)]);
            end
        end
    end
    
    % q_i = min_p q_i^(p) + p
    qi = min([qi, qip(:,p+1) + p],[],2);
    
    % break loop if no lower value for q_i is computable
    % CAUTION: are there cases where this is not a sufficient criterion?
    % ... loop could break if all equations are linear (?)
    if all(qi < p+1)
        break;
    end
    
    % increment power
    p = p + 1;
    Apower = Apower * A;
    
    % safety break condition
    if p == 100
        warning("Check computation for order of abstraction error");
        % fix for prodDesParam, where dot(x4) = 0... thus Inf (probably
        % works also for similar systems)
        if all(qi(~isinf(qi)) < p+1)
            break;
        end
    end
    
end

% write to options struct
options.orders = qi;
options.minorder = min(options.orders);

end

% ------------------------------ END OF CODE ------------------------------
