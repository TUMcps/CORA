function [Rti,Rtp,Rti_y,options] = linReach_adaptive(obj,options,Rstart,Rstart_y)
% linReach_adaptive - computes the reachable set after linearization;
%    automated tuning from [2] applied to algorithm from [1]
%
% Syntax:
%    [Rti,Rtp,options] = linReach_adaptive(obj,options,Rstart)
%
% Inputs:
%    obj - nonlinDASys object
%    options - struct with algorithm settings
%    Rstart - reachable set (state variables)
%    Rstart_y - reachable set (algebraic variables)
%
% Outputs:
%    Rti - reachable set for time interval
%    Rtp - reachable set for time point
%    options - struct with algorithm settings
%
% References:
%    [1] M. Althoff, B. Krogh. "Reachability analysis of nonlinear
%        differential-algebraic systems", IEEE Transactions of
%        Automatic Control, 2013.
%    [2] M. Wetzlinger et al. "Automated parameter tuning for reachability
%        analysis of nonlinear systems", HSCC 2021.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Mark Wetzlinger
% Written:       17-June-2021
% Last update:   31-August-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

error_adm_x = options.error_adm_x_horizon;
error_adm_y = options.error_adm_y_horizon;
% bloating factor
expFactor = 1.1;

% time step size + set abstraction order ----------------------------------
% finite horizon depending on varphi, init at start, control below
lastStep = false; veryfirststep = false; timeStepequalHorizon = false;

% iteration counter run:
% - 0: only for first step in lin, otherwise:
% - 1: compute finitehorizon (use last deltat for varphi estimation),
% - 2: compute optimal deltat

% first step
if options.i == 1
    veryfirststep = true;
    % initial guess for time step size / finite horizon
    options.timeStep = (options.tFinal - options.tStart) * 0.05; % initial guess
    finitehorizon = options.timeStep;
    
    % first step: adaptive tensor order
    options.run = 0;
    options = aux_initStepTensorOrder(obj,options,Rstart,Rstart_y);
    
    % init reduction indices
    options.gredIdx.Rhomti = {};
    options.gredIdx.Rhomtp = {};
    options.gredIdx.Rpar = {};
    
    % init reduction indices (only used for options.tensorOrder = 3)
    options.gredIdx.xyu = {};
    options.gredIdx.Z_error = {};
    options.gredIdx.Z_error_y = {};
    
% last step
elseif options.timeStep > options.tFinal - options.t
    lastStep = true;
    options.timeStep = options.tFinal - options.t;
    finitehorizon = options.timeStep;
    
% non-start/end step
elseif options.i > 1
    % new version: adaptive P-controller
%     finitehorizon = options.finitehorizon(options.i-1) + ...
%         options.kp(options.i-1) * (options.varphi(options.i-1) - options.zetaphi);
    % old version: take last one and estimate so that varphi ~ zetaphi
    finitehorizon = options.finitehorizon(options.i-1) * ...
        (options.decrFactor - options.zetaphi) / ...
        (options.decrFactor - options.varphi(options.i-1));

    % finitehorizon is capped by remaining time
    finitehorizon = min([finitehorizon,options.tFinal - options.t]);

    % init time step size for current step as finitehorizon
    options.timeStep = finitehorizon;
    
    assert(options.timeStep > 0,'Tuning error... report to devs');
end

% check if Rstart is degenerate (each dimension): order may change
zeroWidthDim = sum(abs(generators(zonotope(Rstart))),2) == 0;

options.run = 1;
% -------------------------------------------------------------------------

% log
abscount = 0;

% loop for time step adaptation
while true

    % linearize nonlinDA system
    [obj,linSys,options,linOptions] = linearize(obj,options,Rstart,Rstart_y); 

    if options.run == 1
        % scaling over delta t, required for optimization function
        zetaP = exp(trace(linSys.A*options.timeStep));
    end

    % translate Rstart by linearization point
    Rdelta = Rstart + (-obj.linError.p.x);

    % compute reachable set of linearized system
    % ti: time interval, tp: time point
    linOptions.p = obj.linError.p.x;
    try
        [Rlin,options] = initReach_adaptive(linSys,Rdelta,options,linOptions);
    catch ME % if time step size is too big, then Inf/NaN in expmtie
        if strcmp(ME.identifier,'expmtie:notconverging')
            options.timeStep = options.timeStep * 0.5;
            finitehorizon = options.timeStep;
            continue
        else
            rethrow(ME);
        end
    end
    
    % compute reachable set of the abstracted system including the
    % abstraction error (using the selected evaluation)

    % loop until the actual abstraction error is smaller than
    % the estimated abstraction error
    
    Rlintp = Rlin.tp; Rlinti = Rlin.ti;
    options.Lconverged = false;
    
    perfIndCounter = 1; perfInds = [];
    while true
        % increment counter how often abstraction error computed
        abscount = abscount + 1;

        % convert error to zonotope
        Verror_x = zonotope([0*error_adm_x,diag(error_adm_x)]);
        Verror_y = zonotope([0*error_adm_y,diag(error_adm_y)]);
        % compute Lagrange remainder, see [1, eq.(15)]
        Verror = Verror_x + obj.linError.CF_inv * Verror_y;
        
        % estimate the abstraction error
%         if any(error_adm_x) || any(error_adm_y)
        % compute abstraction error as input solution
        [RallError,options] = errorSolution_adaptive(linSys,options,Verror);
%         disp("step|deltat|run|gens(RallError):" + options.i + "|" + options.timeStep + ...
%              "|" + options.run + "|" + size(generators(RallError),2));
%         else
%             % only used in options.run == 1, but important for adaptive
%             RallError = zonotope(zeros(obj.dim,1));
%         end
        
        try
            % compute overall reachable set including abstraction error
            Rmax = Rlinti + RallError;
            
            % compute abstraction error
            [Verror, errorInt, errorInt_x, errorInt_y, Rti_y, options] = ...
            	abstractionError_adaptive(obj, options, Rmax, Verror_y);
            
        catch ME
            if strcmp(ME.identifier,'reach:setoutofdomain')
                options.Lconverged = false; break
            else
                rethrow(ME);
            end
        end
        
        % compare linearization error with the maximum admissible error
        perfIndCurr = max([errorInt_x ./ error_adm_x; errorInt_y ./ error_adm_y]);
        
        if perfIndCurr <= 1 || ~any(errorInt)
            % error converged... exit
            options.Lconverged = true;
            break
        elseif perfIndCounter > 1
            perfInds(perfIndCounter) = perfIndCurr;
            if perfIndCounter > 2 && perfInds(perfIndCounter) > perfInds(perfIndCounter-1)
                options.Lconverged = false; break
            end
        end
        
        error_adm_x = expFactor * errorInt_x;
        error_adm_y = expFactor * errorInt_y;
        perfIndCounter = perfIndCounter + 1;
    end
    
    
    % if either L (of finitehorizon) was out of domain or did not converge
    if ~options.Lconverged
        options.timeStep = options.timeStep * 0.5;
        finitehorizon = options.timeStep;
        error_adm_x = options.error_adm_x_horizon;
        error_adm_y = options.error_adm_y_horizon;
        continue
    end
    % ... now containment of L ensured
    options.Lconverged = true;
    
    % compute order of abstraction error
    if veryfirststep || ~all(zeroWidthDim == options.zeroWidthDim)
        options.zeroWidthDim = zeroWidthDim;
        options = aux_getPowers(obj,options,linSys,zeroWidthDim,Verror);
    end
    
    % compute set of abstraction errors
    [Rerror,options] = errorSolution_adaptive(linSys,options,Verror);
    
    % measure abstraction error
    abstrerr = sum(abs(generators(Rerror)),2)';
    
    
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
                error_adm_x = options.decrFactor * errorInt_x;
                error_adm_y = options.decrFactor * errorInt_y;
                continue;
            end
            
            % take worst-case current gain (use same dimensions as for lin)
            temp = abstrerr(options.orders == options.minorder) ./ ...
                abstrerr_h(options.orders == options.minorder);
            varphi = max( temp(~isnan(temp)) );
            
            if varphi > options.decrFactor + 10*eps
                throw(CORAerror('CORA:specialError',...
                    'Error in computation of abstraction error'));
            end
            
            % check condition for varphi
            if varphi < options.zetaphi
                % decrease time step size
                finitehorizon = options.timeStep;
                options.timeStep = options.decrFactor * options.timeStep;
                % reset estimates for linearization error
                error_adm_x = options.decrFactor * errorInt_x;
                error_adm_y = options.decrFactor * errorInt_y;
                % svae values from finite horizon
                abstrerr_h = abstrerr;
                Rerror_h = Rerror;
                Rti_h = Rlinti; Rtp_h = Rlintp;
                linx_h = obj.linError.p.x;
                zetaP_h = zetaP;
                continue;
            end
            
            options.varphi(options.i,1) = varphi;
            options.finitehorizon(options.i,1) = finitehorizon;
            % estimate near-optimal delta t from optimization function
            [options.timeStep,~] = aux_optimaldeltat(Rstart,Rerror_h,...
                finitehorizon,varphi,zetaP_h,options);
            % save abstraction error, reset to zero for optimal time step size
            options.error_adm_x_horizon = errorInt_x;
            options.error_adm_y_horizon = errorInt_y;
            error_adm_x = options.error_adm_x_Deltatopt;
            error_adm_y = options.error_adm_y_Deltatopt;
            
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
            % save abstraction error
            options.error_adm_x_horizon = errorInt_x;
            options.error_adm_y_horizon = errorInt_y;
            % set abstraction error to values for previous optimal delta t
            error_adm_x = options.error_adm_x_Deltatopt;
            error_adm_y = options.error_adm_y_Deltatopt;
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
                options.timeStep,Rerror_h,Rerror,options.decrFactor,options.minorder);
        end
        % save finite horizon and abstraction error for next step
        options.finitehorizon(options.i,1) = finitehorizon;
        options.error_adm_x_Deltatopt = errorInt_x;
        options.error_adm_y_Deltatopt = errorInt_y;
        
        % predict tensor order of next step (unless fixed)
        if ~options.fixedTensorOrder
            if options.i == 1
                % save time step and abstrerr
                options.kappa_deltat = options.timeStep;
                options.kappa_abstrerr = vecnorm(abstrerr);
            
                % compute finitehorizon of next step acc. to extrapolation
%                 temp = options.finitehorizon(options.i) * ...
%                     (options.decrFactor - options.zetaphi) / ...
%                     (options.decrFactor - options.varphi(options.i));
                % compute gain kp
%                 options.kp(options.i,1) = ...
%                     (temp - options.finitehorizon(options.i)) / ...
%                     (options.varphi(options.i) - options.zetaphi);
            
            elseif ~lastStep
                % compensate for the difference between the current time
                % step size and stored time step size by varphi
                options = aux_nextStepTensorOrder(obj,options,...
                    vecnorm(abstrerr),linSys,Rmax,Verror_y);
            
                options = aux_nextStepTensorOrder_new(obj,options,...
                    abstrerr,linSys,Rmax,Verror_y,Rlintp);
            end
            
            % update gain (how much gain kp would have been necessary to
            % achieve varphi = zetaphi in this step?)
%             options.kp(options.i,1) = ...
%                 (options.finitehorizon(options.i) - options.finitehorizon(options.i-1)) / ...
%                 (options.varphi(options.i) - options.zetaphi);
            
        end
        
        % second run finished -> exit
        break
    end
    
    
    % update counter
    options.run = options.run + 1;
    
    if options.timeStep == finitehorizon
        % required Rerror and Rlin already there, but compute
        % timeStep = decrFactor*finitehorizon for usable varphi
        % (otherwise two points required for varphi are the same point)
        timeStepequalHorizon = true;
        options.timeStep = options.timeStep * options.decrFactor;
    end 
    
end


% use correct Rerror in case optimal time step size is equal to finitehorizon
if timeStepequalHorizon
    % reset time step size for time increment in loop outside
    options.timeStep = finitehorizon;
    % choose sets computed by Delta (last run: finitehorizon*decrFactor)
    Rti = Rti_h + linx_h; Rtp = Rtp_h + linx_h;
    Rerror = Rerror_h;
else
    % translate reachable sets by linearization point
    Rti = Rlinti + obj.linError.p.x;
    Rtp = Rlintp + obj.linError.p.x;
end

% add the abstraction error to the reachable sets
Rti = Rti + Rerror;
Rtp = Rtp + Rerror;

% save counter of how often Lagrange remainder was computed
options.abscount(options.i,1) = abscount;

end


% Auxiliary functions -----------------------------------------------------

% adaptation of tensorOrder: init step
function options = aux_initStepTensorOrder(obj,options,Rstartset,R_y)
% computes the lower bound for the Lagrange remainder
% ... that is, L and Rerr for the start set of the current step

options.fixedTensorOrder = false;
if isfield(options,'tensorOrder') && any(options.tensorOrder == [2,3])
    options.fixedTensorOrder = true;
    return;
end

% linearization point u* of the input is the center of the current input set
obj.linError.p.u = center(options.U) + options.uTrans;

% linearization points in differential set (x*) and algebraic set (y*)
obj.linError.p.x = center(Rstartset);
obj.linError.p.y = consistentInitialState(obj, obj.linError.p.x, center(R_y), obj.linError.p.u);

% substitute p into the system equation to obtain the constant input
f0_dyn = obj.dynFile(obj.linError.p.x, obj.linError.p.y, obj.linError.p.u); % f(z*) in [1]
f0_con = obj.conFile(obj.linError.p.x, obj.linError.p.y, obj.linError.p.u); % g(z*) in [1]

% evaluate Jacobians
[A,B,C,D,E,F] = obj.jacobian(obj.linError.p.x, obj.linError.p.y, obj.linError.p.u);

% matrices of the linearized system (see [1, eq.(14)])
obj.linError.D = D;
obj.linError.E = E;
obj.linError.F_inv = pinv(F);
obj.linError.CF_inv = C*obj.linError.F_inv;
obj.linError.f0 = f0_dyn - obj.linError.CF_inv*f0_con;
obj.linError.f0_con = f0_con;

% shift start set by linearization point
Rdelta = Rstartset + (-obj.linError.p.x);
% estimate at the beginning of the loop also all-zero
Verror_y = zonotope(zeros(obj.nrOfConstraints,1));

% compute L for both kappas
options.tensorOrder = 2;
[~, L_2, L2_x, L2_y] = abstractionError_adaptive(obj, options, Rdelta, Verror_y);
options.tensorOrder = 3;
[~, L_3, L3_x, L3_y] = abstractionError_adaptive(obj, options, Rdelta, Verror_y);

% remove linear dimensions for comparison
L_2lin = L_2(L_2 ~= 0);
L_3lin = L_3(L_3 ~= 0);
% Lxy_2lin = [L2_x; L2_y];
% Lxy_2lin = Lxy_2lin(Lxy_2lin ~= 0);
% Lxy_3lin = [L3_x; L3_y];
% Lxy_3lin = Lxy_3lin(Lxy_3lin ~= 0);

% compare Lagrange remainder to decide which tensorOrder to use
if all( L_3lin ./ L_2lin > options.zetaK )
% if all( Lxy_3lin ./ Lxy_2lin > options.zetaK )
    options.tensorOrder = 2;
else
    options.tensorOrder = 3;
end

end

% adaptation of tensorOrder: every step for next step
function options = aux_nextStepTensorOrder(obj,options,abstrerr,linSys,Rmax,Verror_y)

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
    % only check other tensorOrder if abstraction error has grown
    if abs(abstrerr_est / abstrerr_prev - 1) > 1 - options.zetaK
        abstrerr_2 = abstrerr;
        % compute abstraction error for other tensor order
        options.tensorOrder = 3;
    else
        return;
    end
    
    % compute set of abstraction errors
    Verror = abstractionError_adaptive(obj, options, Rmax, Verror_y);
    Rerror = errorSolution_adaptive(linSys,options,Verror);
    % simplify to scalar value
    err_3 = sum(abs(generators(Rerror)),2);
    abstrerr_3 = vecnorm(err_3);
    
    % compare results
    if ~all( abstrerr_3 ./ abstrerr_2 > options.zetaK )
        options.kappa_deltat = options.timeStep;
        options.kappa_abstrerr = abstrerr_3;
        options.error_adm_Deltatopt = err_3;
        options.error_adm_horizon = zeros(obj.dim,1);
    else
        options.tensorOrder = 2;
    end
    
elseif options.tensorOrder == 3
    % only check other tensorOrder if abstraction error has shrunk
    if abs(abstrerr_est / abstrerr_prev - 1) > 1 - options.zetaK
        abstrerr_3 = abstrerr;
        % compute other tensor order and check again
        options.tensorOrder = 2;
    else
        return;
    end
    
    % compute set of abstraction errors
    Verror = abstractionError_adaptive(obj,options,Rmax,Verror_y);
    Rerror = errorSolution_adaptive(linSys,options,Verror);
    % simplify to scalar value
    err_2 = sum(abs(generators(Rerror)),2);
    abstrerr_2 = vecnorm(err_2);
    
    % compare results
    if all( abstrerr_3 ./ abstrerr_2 > options.zetaK )
        options.kappa_deltat = options.timeStep;
        options.kappa_abstrerr = abstrerr_2;
        options.error_adm_Deltatopt = err_2;
        options.error_adm_horizon = zeros(obj.dim,1);
    else
        options.tensorOrder = 3;
    end
end

end

% adaptation of tensorOrder: every step for next step
function options = aux_nextStepTensorOrder_new(obj,options,abstrerr,linSys,Rmax,Verror_y,Rlintp)

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
    % only check other tensorOrder if abstraction error has grown
    if abs(abstrerr_est / abstrerr_prev - 1) > 1 - options.zetaK
        abstrerr_2 = abstrerr';
        % compute abstraction error for other tensor order
        options.tensorOrder = 3;
    else
        return;
    end
    
    % compute set of abstraction errors
    Verror = abstractionError_adaptive(obj, options, Rmax, Verror_y);
    Rerror = errorSolution_adaptive(linSys,options,Verror);
    % compute diameter
    abstrerr_3 = sum(abs(generators(Rerror)),2);
    % compute cost values
    diamRtp = sum(abs(generators(Rlintp)),2);
    remSteps = (options.tFinal - options.t) / timeStep_prev;
    estval_2 = (1 + abstrerr_2 ./ diamRtp) .^ remSteps;
    estval_3 = (1 + abstrerr_3 ./ diamRtp) .^ remSteps;
    
    % compare results
    if ~all( estval_3 ./ estval_2 > options.zetaK )
        options.kappa_deltat = options.timeStep;
        options.kappa_abstrerr = vecnorm(abstrerr_3);
        options.error_adm_Deltatopt = abstrerr_3;
        options.error_adm_horizon = zeros(obj.dim,1);
    else
        options.tensorOrder = 2;
    end
    
elseif options.tensorOrder == 3
    % only check other tensorOrder if abstraction error has shrunk
    if abs(abstrerr_est / abstrerr_prev - 1) > 1 - options.zetaK
        abstrerr_3 = abstrerr';
        % compute other tensor order and check again
        options.tensorOrder = 2;
    else
        return;
    end
    
    % compute set of abstraction errors
    Verror = abstractionError_adaptive(obj,options,Rmax,Verror_y);
    Rerror = errorSolution_adaptive(linSys,options,Verror);
    % simplify to scalar value
    abstrerr_2 = sum(abs(generators(Rerror)),2);
    % compute cost values
    diamRtp = sum(abs(generators(Rlintp)),2);
    remSteps = (options.tFinal - options.t) / timeStep_prev;
    estval_2 = (1 + abstrerr_2 ./ diamRtp) .^ remSteps;
    estval_3 = (1 + abstrerr_3 ./ diamRtp) .^ remSteps;
    
    % compare results
    if all( estval_3 ./ estval_2 > options.zetaK )
        options.kappa_deltat = options.timeStep;
        options.kappa_abstrerr = vecnorm(abstrerr_2);
        options.error_adm_Deltatopt = abstrerr_2;
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
dHused = 0.5;
zetaZ = opt.redFactor * dHused; % better representation of actual reduction

% kprimemax = log(deltatmin/deltat) / log(mu); % if deltatmin computed
kprimemax = ceil(-log(100) / log(mu)); % max 100 steps instead of 1 (10000 total steps max)
kprime = 0:round(kprimemax); % just to catch rounding errors
k = mu .^ (-kprime);
deltats = deltat * mu.^kprime; % consecutively scaled down by mu
% number of "full" steps, i.e., steps where deltats(i) fully used (all but last)
floork = floor(k);

% set size estimates
rR = vecnorm(sum(abs(generators(Rt)),2),2);
rerr1 = vecnorm(sum(abs(generators(Rerr)),2),2);

% model varphi by linear interpolation over the gain
% ...bounds: 1 (at deltatmin) to measured gain (at deltat)
varphimax = mu; varphiDelta = (varphimax - varphimin);
varphi = (varphimin + (deltats(1) - deltats)/deltats(1) * varphiDelta) / mu;
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
function varphi = aux_varphiest(horizon,deltat,Rerr_h,Rerr_deltat,decrFactor,minorder)

% compute estimate for set sizes
if isa(Rerr_h,'zonotope') && isa(Rerr_deltat,'zonotope')
    rerr1 = vecnorm(sum(abs(generators(Rerr_h)),2),2);
    rerrk = vecnorm(sum(abs(generators(Rerr_deltat)),2),2);
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
residual_up = varphitotal - rhs;

% lower bound: varphi = 0
varphi_low = 0;
% compute accumulated scaling
varphitotal = aux_varphitotal(varphi_low,varphi_lim,deltat,horizon,decrFactor);
% compute residual of prod of varphi (assuming varphi) and rhs
residual_low = varphitotal - rhs;

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
        residual_low = residual(cnt);
        varphi_low = varphi(cnt);
    else
        residual_up = residual(cnt);
        varphi_up = varphi(cnt);
    end
    
    % break condition
    if cnt == 10000
        throw(CORAerror('CORA:notConverged','Internal estimation'));
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
% generally, only called once after the very first converged L!

% hard-coded shortcut
if options.i == 1 && isfield(options,'orders')
    options.minorder = min(options.orders);
    return;
end

% propagation matrix
A = linSys.A;
n = obj.dim;
    
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
    if isfield(options,'isHessianConst') && options.isHessianConst
        % easier case: Hessian is constant
        Hess = options.hessianConst;
    else
        if strcmp(options.alg,'poly')
            Hess = obj.hessian(obj.linError.p.x,obj.linError.p.y,obj.linError.p.u);
        else
            % Hessian requires evaluation on intervals for x and u...
            % try eval on linearization point, set handle to correct file
            obj = setHessian(obj,'standard');
            Hess = obj.hessian(obj.linError.p.x,obj.linError.p.y,obj.linError.p.u);
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
    
    % safety condition
    if p == 100
        throw(CORAerror('CORA:notConverged','Order'));
    end
    
end

% write to options struct
options.orders = qi;
options.minorder = min(options.orders);

end

% ------------------------------ END OF CODE ------------------------------
