classdef linErrorBound < handle
% linErrorBound - helper class for dealing with error bounds in the
%    adaptive reachability algorithm for linear continuous-time systems
% note: do not use outside of built-in reachability algorithms
%
% Syntax:
%    errs = linErrorBound(emax,tFinal)
%
% Inputs:
%    emax - error bound
%    tFinal - time horizon
%
% Outputs:
%    errs - generated linErrorBound object
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
% See also: linearSys/private/priv_reach_adaptive

% Authors:       Mark Wetzlinger
% Written:       06-November-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = public, GetAccess = public)  % Access = private
    
    % --- errors ---
    emax = [];              % full error margin
    tFinal = [];            % time horizon
    timeSteps = {};         % time step sizes

    step_acc = {};          % committed accumulating error (k-th step)
    step_red = [];          % committed reduction error (k-th step)

    seq_nonacc = {};        % committed non-accumulating error (each step)
    cum_acc = [];           % committed accumulating error (cumulative)
    cum_red = [];           % committed reduction error (cumulative)

    bound_acc = [];         % bound for accumulating error
    bound_red = [];         % bound for reduction error
    bound_red_max = [];     % maximum value of reduction error (at tFinal)
                            % (note: bound_acc_max = emax - bound_red_max)
    bound_rem = [];         % bound for remaining error (after acc/red)
    bound_remacc = [];      % bound for remaining accumulating error

    bound_acc_ok = {};      % fulfillment of accumulating error bounds
    bound_nonacc_ok = {};   % fulfillment of non-accumulating error bounds

    idv_PUtkplus1 = {};     % individual error: particular solution (accumulating)
    idv_F = {};             % individual error: state curvature (non-accumulating)
    idv_G = {};             % individual error: input curvature (non-accumulating)
    idv_linComb = {};       % individual error: linear combination (non-accumulating)
    idv_PUtauk = {};        % individual error: time dependency of particular solution (non-accumulating)

    useApproxFun = [];      % flag whether to use approximation functions (true)
                            % or bisection (false) to estimate time step size

    % --- coefficients for approximation functions ---
    % order of the individual terms over Delta t
    % - eps_linComb:    linear
    % - eps_F:          quadratic
    % - eps_G:          quadratic
    % - eps_PU:         quadratic
    
    % coefficients for approximation functions
    % - linear:         eps(Delta t) =                 b * Delta t
    % - quadratic:      eps(Delta t) = a * Delta t^2 + b * Delta t
    % ... no offset since all errors are 0 for Delta t = 0
    coeff_PUtkplus1_a = [];
    coeff_PUtkplus1_b = [];
    coeff_F_a = [];
    coeff_F_b = [];
    coeff_G_a = [];
    coeff_G_b = [];
    coeff_linComb_b = [];
    coeff_PUtauk_b = [];

    % for bisection
    bisect_lb_timeStep_acc = [];
    bisect_lb_acc = [];
    bisect_lb_accok = [];
    bisect_lb_acc_perc = [];
    bisect_lb_timeStep_nonacc = [];
    bisect_lb_nonacc = [];
    bisect_lb_nonaccok = [];
    bisect_lb_nonacc_perc = [];

    bisect_ub_timeStep_acc = [];
    bisect_ub_acc = [];
    bisect_ub_accok = [];
    bisect_ub_acc_perc = [];
    bisect_ub_timeStep_nonacc = [];
    bisect_ub_nonacc = [];
    bisect_ub_nonaccok = [];
    bisect_ub_nonacc_perc = [];

end

methods
    % constructor
    function errs = linErrorBound(emax,tFinal)
        % check number of input arguments
        narginchk(2,2);
        
        % check input arguments
        if ~isnumeric(emax) || ~isscalar(emax) || ~isfinite(emax) || emax <= 0
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Error margin must be a scalar, positive, real number.'));
        end
        if ~isnumeric(tFinal) || ~isscalar(tFinal) || ~isfinite(tFinal) || tFinal < 0
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Time horizon must be a scalar, non-negative, real number.'));
        end

        errs.emax = emax;
        errs.tFinal = tFinal;
        % ...all other values initialized by empty (see properties)
    end

    function res = checkErrors(errs,Rout_error,Rout_tp_error)
        % check function which returns true/false depending on fulfillment
        % of error bounds over time

        % assume checks to be ok
        res = true;
        internalcheck = true;

        % read out vector of time steps
        tVec = cell2mat(errs.timeSteps);
        
        % read out some errors
        nonacc_step = cell2mat(errs.seq_nonacc);
        acc_step = cell2mat(errs.step_acc);
        nonacc_bound = cell2mat(errs.bound_rem);
        acc_bound = cell2mat(errs.bound_acc);
        red_bound = cell2mat(errs.bound_red);

        % check which steps are computed completely
        fullstep = ~isnan(nonacc_step);

        % quick fix for time-point error (empty in verify-call)
        if isempty(Rout_tp_error)
            Rout_tp_error = zeros(numel(tVec),1);
        end
        
        % cumulative sum to yield bound for accumulating error
        % note: subtract error committed in the step!
        acctotal_bound = errs.cum_acc + acc_bound - acc_step;
        % redtotal_bound = errs.cum_red + red_bound - ered;
        
        % compute total committed error
        total_error = nonacc_step + errs.cum_acc - acc_step + errs.cum_red;
        total_error_tp = errs.cum_acc + errs.cum_red;
        
        % --- checks ---
        % all errors and bounds need to be larger or equal to zero
        if any([Rout_error(fullstep); Rout_tp_error(fullstep); nonacc_step(fullstep); ...
                acc_step; errs.step_red; nonacc_bound; acc_bound; red_bound; acctotal_bound] < 0)
            res = false;
        end
        
        % full errors (time-interval and time-point)
        % (1) need to be below maximum error
        if any(Rout_error(fullstep) > errs.emax) || any(Rout_tp_error(fullstep) > errs.emax)
            res = false;
        end
        % (2) re-computation of full errors has to match ongoing computation
        %     note: use a tolerance here...
        if any(abs(total_error(fullstep) - Rout_error(fullstep)) > 1e-9) || ...
                (any(Rout_tp_error) && any(abs(total_error_tp(fullstep) - Rout_tp_error(fullstep)) > 1e-9))
            internalcheck = false;
        end
        
        % non-accumulating errors (linComb, F, G)
        % (1) need to be below maximum error
        % (2) need to be below remaining error (after subtraction of acc error
        %     until now and reduction error bound)
        if any(nonacc_step > errs.emax)
            res = false;
        end
        if any(nonacc_step > nonacc_bound)
            internalcheck = false;
        end
        
        % accumulating error (PU)
        % (1) single step error needs to be below accumulated error
        % (2) single step error needs to be below corresponding bound
        % (3) accumulated error needs to be below maximum error
        % (4) accumulated error needs to be below final admissible value
        if any(errs.cum_acc > errs.emax) || any(errs.cum_acc > errs.emax - errs.bound_red_max)
            res = false;
        end
        if any(acc_step > errs.cum_acc) || any(acc_step > acc_bound) 
            internalcheck = false;
        end
        % (5) accumulated error needs to be below linearly increasing bound
        %     (note: only in the current version!)
        % (6) accumulated error bound needs to be below linearly increasing bound
        % compute bound first...
        finalValue = errs.emax - errs.bound_red_max;
        cumsumtVec_tp = cumsum(tVec);
        linBound = finalValue * cumsumtVec_tp / cumsumtVec_tp(end);
        if any(errs.cum_acc > linBound) || any(linBound - acctotal_bound < -10*eps)
            internalcheck = false;
        end
        
        % reduction error (PU)
        % (1) single step error needs to be below accumulated error
        % (2) single step error needs to be below corresponding bound
        % (3) accumulated error needs to be below maximum error
        % (4) accumulated error needs to be below final admissible value
        if any(errs.cum_red > errs.emax)
            res = false;
        end
        if any(errs.step_red > errs.cum_red) || any(errs.step_red > red_bound) || ...
                 any(errs.cum_red > errs.bound_red_max)
            internalcheck = false;
        end
        
        % print warning if internal checks failed
        if ~internalcheck
            CORAwarning("CORA:contDynamics","Internal checks failed");
        end

    end

    function errs = nextBounds(errs,timeStep,t,k,k_iter)
        % compute error bounds given a time step size
        errs.timeSteps{k,1}(k_iter,1) = timeStep;

        % compute combined error margin for eacc and enonacc this step:
        % given by the maximum error minus the final value for the
        % reduction error and the sum of all accumulating errors until
        % the current step
        if k == 1
            errs.bound_remacc(k,1) = errs.emax - errs.bound_red_max;
        else
            errs.bound_remacc(k,1) = errs.emax - errs.bound_red_max - errs.cum_acc(k-1);
        end

        % the accumulating error bound can be at most the value of the
        % combined error margin scaled by the fraction the time step size
        % covers of the remaining time
        errs.bound_acc{k,1}(k_iter,1) = errs.bound_remacc(k,1) * timeStep / (errs.tFinal - t);

        % read out value for admissible reduction error from curve:
        % use the corresponding linearly interpolated error value as the 
        % admissible error (also, subtract the reduction error until now, 
        % since the adaptive-reduce function only considers *local* bounds)
        bound_red_tk = errs.bound_red_max * (t+timeStep) / errs.tFinal;
        if k == 1
            errs.bound_red{k,1}(k_iter,1) = bound_red_tk;
            errs.bound_rem{k,1}(k_iter,1) = errs.emax - bound_red_tk;
        else
            errs.bound_red{k,1}(k_iter,1) = bound_red_tk - errs.cum_red(k-1);
            errs.bound_rem{k,1}(k_iter,1) = errs.emax - bound_red_tk - errs.cum_acc(k-1);
        end

    end

    function accumulateErrors(errs,k,k_iter)
        % compute accumulating errors (acc/red), also save non-accumulating
        % errors in sequence vector

        if k == 1
            errs.cum_acc = errs.step_acc{k}(k_iter);
            errs.cum_red = errs.step_red;
        else
            errs.cum_acc(end+1,1) = errs.cum_acc(end) + errs.step_acc{k}(k_iter);
            errs.cum_red(end+1,1) = errs.cum_red(end) + errs.step_red(k);
        end

    end

    function removeRedundantValues(errs,k,k_iter)
        % in each step k, we save the errors for all iterations k_iter
        % however, after a specific time step size is chosen, we only keep
        % the error which has actually been committed
        errs.seq_nonacc{k,1} = errs.seq_nonacc{k}(k_iter);

        errs.idv_PUtkplus1{k,1} = errs.idv_PUtkplus1{k}(k_iter);
        errs.idv_F{k,1} = errs.idv_F{k}(k_iter);
        errs.idv_G{k,1} = errs.idv_G{k}(k_iter);
        errs.idv_linComb{k,1} = errs.idv_linComb{k}(k_iter);
        errs.idv_PUtauk{k,1} = errs.idv_PUtauk{k}(k_iter);

        errs.step_acc{k,1} = errs.step_acc{k}(k_iter);

        errs.timeSteps{k,1} = errs.timeSteps{k}(k_iter);

        errs.bound_rem{k,1} = errs.bound_rem{k}(k_iter);
        errs.bound_acc{k,1} = errs.bound_acc{k}(k_iter);
        errs.bound_red{k,1} = errs.bound_red{k}(k_iter);

    end

    function [Rcont_error,Rcont_tp_error] = fullErrors(errs,k)
        % compute the errors of the time-point and time-interval reachable
        % set; since the non-accumulating error contains the accumulating
        % error of the last step, we have two cases:
        if k == 1
            Rcont_error = errs.seq_nonacc{k} + errs.cum_red(k);
        else
            Rcont_error = errs.seq_nonacc{k} + errs.cum_acc(k-1) + errs.cum_red(k);
        end
        Rcont_tp_error = errs.cum_acc(k) + errs.cum_red(k);
    end

    function eps_linComb = compute_eps_linComb(errs,eAdt,startset)
        % compute error of linear combination, see [1, Prop. 9]
        n = size(eAdt,1);
        G_minus = (eAdt - eye(n)) * generators(startset);
        eps_linComb = sqrt( size(G_minus,2) ) * norm(G_minus,2);
    end

    function [eps_F,box_Fstartset_c,box_Fstartset_G] = compute_eps_F(errs,F,startset)
        % compute curvature error (state), see [1, Prop. 1]:
        %    eps_F = 2 * errOp(F * startset)

        if isa(startset,'zonotope')
            % small speed up for zonotopes
            Fcenter = center(F); Frad = rad(F);
            box_Fstartset_c = Fcenter * startset.c;
            Fstartset_G = [Fcenter * startset.G, ...
                diag(Frad * sum(abs([startset.c, startset.G]),2))];
            box_Fstartset_G = sum(abs(Fstartset_G),2);
            eps_F = 2 * vecnorm(box_Fstartset_G + abs(box_Fstartset_c));
        else
            % standard computation
            errorset = F*startset;
            eps_F = 2 * priv_errOp(errorset);
            if nargout > 1
                Z_errorset = zonotope(errorset);
                box_Fstartset_c = Z_errorset.c;
                box_Fstartset_G = Z_errorset.G;
            end
        end
    end

    function [eps_G,Gu_c,Gu_G] = compute_eps_G(errs,G,u)
        % compute curvature error (input), see [1, Prop. 1]:
        %    eps_G = 2 * errOp(G*u)
        Gu_c = center(G) * u;
        Gu_G = rad(G) * abs(u);
        eps_G = 2 * vecnorm(Gu_G + abs(Gu_c));
    end

    function eps_PUtauk = compute_eps_PUtauk(errs,eAt,PU)
        % compute error due to loss of time dependency when adding the
        % particular solution to the homogeneous solution, see [1, Prop. 3]
        eps_PUtauk = priv_errOp(eAt*PU);
    end

    function eps_PUtkplus1 = compute_eps_PUtkplus1(errs,eAt,U,Asum,Asum_U,E,timeStep)
        % compute error in particular solution, see [1, Prop. 2]
        eps_PUtkplus1 = priv_errOp(eAt*(Asum*U + E*timeStep*U)) ...
            + priv_errOp(eAt*(Asum_U+E*timeStep*U));
    end

    function computeErrorBoundReduction(errs,A,G_U)
        % implementation of heuristics from [1, Sec. IV.D.2)]:
        % determine a curve
        %   x-axis (time):  0 (tStart)        -> tFinal
        %   y-axis (error): 0 (initial error) -> bound_red_max
        % for the reduction error to yield the smallest zonotope order at
        % the end of the time horizon; note that we only require to reduce
        % the particular solution due to the uncertainty in the input (U);
        % the result is a curve from which we can read the admissible
        % reduction error and thus also the remaining error margin for the
        % other errors for any point in time

        % no generators
        if ~any(G_U,'all')
            errs.bound_red_max = 0;
            return
        end

        % heuristics for near-optimal allocation
        stepsforerrorbound = 100;
        timeStep = errs.tFinal / stepsforerrorbound;
        n = size(A,1);
        e_At = eye(n);
        e_ADeltatk = expm(A * timeStep);
        
        % 1. compute auxiliary sets V_k and errors e(V_k)
        V = cell(stepsforerrorbound,1);
        errV = zeros(stepsforerrorbound,1);
        
        DeltatkU = timeStep * G_U;
        for i=1:stepsforerrorbound
            % compute auxiliary set
            V{i} = e_At * DeltatkU;
            % propagate exponential matrix for next V
            e_At = e_At * e_ADeltatk;
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
        bestfinalorder = stepsforerrorbound+1;
        min_idx = 1;
        
        % loop over the mesh of emax to find a near-optimal result
        for i=1:meshsize-1
            % find chi* = number of V_k that can be reduced for current red error
            idx = find(errVsort_cumsum < errs.emax * emax_percentage4ered(i),1,'last');
            % fraction emax_percentage4ered is allocated to the reduction error,
            % results in factor N of total number of steps
            N = 1/(1-emax_percentage4ered(i));
            % compute resulting zonotope order
            if isempty(idx)
                finalorder(i) = N*(stepsforerrorbound+1);
            else
                finalorder(i) = N*(stepsforerrorbound+1 - idx);
            end
            % best reduction error allocation is the one which results in the
            % lowest zonotope order, also save the corresponding idx
            if finalorder(i) < bestfinalorder
                bestfinalorder = finalorder(i);
                min_idx = idx;
            end
        end
        
        % simpler method
        errs.bound_red_max = errs.emax * emax_percentage4ered(min_idx);
    end

    function updateCoefficientsApproxFun(errs,k,k_iter,fullcomp)
    % compute the coefficients for the approximation functions which model
    % the behavior of the individual error over Delta t

    % only if approximation function should be used
    if ~errs.useApproxFun
        return
    end

    if k_iter == 1
        % initialize coefficients (guessing that b = 0)
        
        if fullcomp
            % linComb and time-interval error of PU are linear approximation
            % functions
            errs.coeff_linComb_b = errs.idv_linComb{k}(k_iter) / errs.timeSteps{k}(k_iter);
            errs.coeff_PUtauk_b = errs.idv_PUtauk{k}(k_iter) / errs.timeSteps{k}(k_iter);
        
            % all others are quadratic approximation functions
            % here, we guess that b = 0
            errs.coeff_F_a = errs.idv_F{k}(k_iter) / errs.timeSteps{k}(k_iter)^2;
            errs.coeff_F_b = 0;
            errs.coeff_G_a = errs.idv_G{k}(k_iter) / errs.timeSteps{k}(k_iter)^2;
            errs.coeff_G_b = 0;
        end
        
        errs.coeff_PUtkplus1_a = errs.idv_PUtkplus1{k}(k_iter) / errs.timeSteps{k}(k_iter)^2;
        errs.coeff_PUtkplus1_b = 0;

    else
        % update coefficients
        
        if fullcomp
            % linear approximation function -> take newest value
            errs.coeff_linComb_b = errs.idv_linComb{k}(k_iter) / errs.timeSteps{k}(k_iter);
            errs.coeff_PUtauk_b = errs.idv_PUtauk{k}(k_iter) / errs.timeSteps{k}(k_iter);
        end
        
        % quadratic approximation functions -> take most recent two values
        % a * last_timeStep^2 + b * last_timeStep = lastEps
        % a * timeStep^2      + b * timeStep      = eps
        % -> Deltatmat * [a;b] = eps ... solve using '\'-operator
        Deltatmat = [errs.timeSteps{k}(k_iter-1)^2, errs.timeSteps{k}(k_iter-1); ...
                     errs.timeSteps{k}(k_iter)^2, errs.timeSteps{k}(k_iter)];
    
        % sanity check of singularity of Deltatmat
        if abs(1/cond(Deltatmat)) < eps
            throw(CORAerror('CORA:notConverged','Estimation of time step size'));
        end
    
        % update coefficient of approximation function for eps_PU
        coeffs_PUtkplus1 = Deltatmat \ [errs.idv_PUtkplus1{k}(k_iter-1); errs.idv_PUtkplus1{k}(k_iter)];
        errs.coeff_PUtkplus1_a = coeffs_PUtkplus1(1);
        errs.coeff_PUtkplus1_b = coeffs_PUtkplus1(2);
        
        if fullcomp
            % update coefficient of approximation function for eps_F
            coeffs_F = Deltatmat \ [errs.idv_F{k}(k_iter-1); errs.idv_F{k}(k_iter)];
            errs.coeff_F_a = coeffs_F(1);
            errs.coeff_F_b = coeffs_F(2);
            
            if errs.coeff_F_b < 0
                % then a certain region is smaller than 0 which cannot ever happen
                % -> use only current value and set errs.coeff_F.b to 0
                errs.coeff_F_b = 0;
                errs.coeff_F_a = errs.idv_F{k}(k_iter) / errs.timeSteps{k}(k_iter)^2;
            end
        
            % update coefficient of approximation function for eps_G
            coeffs_G = Deltatmat \ [errs.idv_G{k}(k_iter-1); errs.idv_G{k}(k_iter)];
            errs.coeff_G_a = coeffs_G(1);
            errs.coeff_G_b = coeffs_G(2);
        end

    end

    end

    function updateBisection(errs,k,k_iter,isU,timeStep)

        if k_iter == 1

            % initialize lower bound with 0
            errs.bisect_lb_timeStep_acc = 0;
            errs.bisect_lb_acc = 0;
            errs.bisect_lb_accok = true;
            errs.bisect_lb_acc_perc = 0;
            
            errs.bisect_lb_timeStep_nonacc = 0;
            errs.bisect_lb_nonacc = 0;
            errs.bisect_lb_nonaccok = true;
            errs.bisect_lb_nonacc_perc = 0;
            
            % upper bound is time step size
            errs.bisect_ub_timeStep_acc = timeStep;
            errs.bisect_ub_acc = errs.idv_PUtkplus1{k}(k_iter);
            errs.bisect_ub_accok = errs.bound_acc_ok{k}(k_iter);
            errs.bisect_ub_acc_perc = errs.step_acc{k}(k_iter) / errs.bound_acc{k}(k_iter);
            
            errs.bisect_ub_timeStep_nonacc = timeStep;
            errs.bisect_ub_nonacc = errs.seq_nonacc{k}(k_iter);
            errs.bisect_ub_nonaccok = errs.bound_nonacc_ok{k}(k_iter);
            errs.bisect_ub_nonacc_perc = errs.seq_nonacc{k}(k_iter) / errs.bound_rem{k}(k_iter);

        else

            % determine whether new time step size is new lb or new ub
            if isU
                if ~errs.bound_acc_ok{k}(k_iter) || timeStep > errs.bisect_ub_timeStep_acc
                    % current time step size is always new ub if errcheck not ok
                    
                    if timeStep > errs.bisect_ub_timeStep_acc
                        % current ub becomes lb
                        errs.bisect_lb_timeStep_acc = errs.bisect_ub_timeStep_acc;
                        errs.bisect_lb_acc = errs.bisect_ub_acc;
                        errs.bisect_lb_accok = errs.bisect_ub_accok;
                        errs.bisect_lb_acc_perc = errs.bisect_ub_acc_perc;
                    end
                    % assign new ub
                    errs.bisect_ub_timeStep_acc = timeStep;
                    errs.bisect_ub_acc = errs.idv_PUtkplus1{k}(k_iter);
                    errs.bisect_ub_accok = errs.bound_acc_ok{k}(k_iter);
                    errs.bisect_ub_acc_perc = errs.step_acc{k}(k_iter) / errs.bound_acc{k}(k_iter);
                else % smaller than before and errcheck ok -> new lb
                    errs.bisect_lb_timeStep_acc = timeStep;
                    errs.bisect_lb_acc = errs.idv_PUtkplus1{k}(k_iter);
                    errs.bisect_lb_accok = errs.bound_acc_ok{k}(k_iter);
                    errs.bisect_lb_acc_perc = errs.step_acc{k}(k_iter) / errs.bound_acc{k}(k_iter);
                end
            end
            
            if ~errs.bound_nonacc_ok{k}(k_iter) || timeStep > errs.bisect_ub_timeStep_nonacc
                % current time step size is always new ub if errcheck not ok
                
                if timeStep > errs.bisect_ub_timeStep_nonacc
                    % current ub becomes lb
                    errs.bisect_lb_timeStep_nonacc = errs.bisect_ub_timeStep_nonacc;
                    errs.bisect_lb_nonacc = errs.bisect_ub_nonacc;
                    errs.bisect_lb_nonaccok = errs.bisect_ub_nonaccok;
                    errs.bisect_lb_nonacc_perc = errs.bisect_ub_nonacc_perc;
                end
                errs.bisect_ub_timeStep_nonacc = timeStep;
                errs.bisect_ub_nonacc = errs.seq_nonacc{k}(k_iter);
                errs.bisect_ub_nonaccok = errs.bound_nonacc_ok{k}(k_iter);
                errs.bisect_ub_nonacc_perc = errs.seq_nonacc{k}(k_iter) / errs.bound_rem{k}(k_iter);
            else % smaller than before and errcheck ok -> new lb
                errs.bisect_lb_timeStep_nonacc = timeStep;
                errs.bisect_lb_nonacc = errs.seq_nonacc{k}(k_iter);
                errs.bisect_lb_nonaccok = errs.bound_nonacc_ok{k}(k_iter);
                errs.bisect_lb_nonacc_perc = errs.seq_nonacc{k}(k_iter) / errs.bound_rem{k}(k_iter);
            end

        end

    end

    function timeStep = estimateTimeStepSize(errs,t,k,k_iter,fullcomp,timeStep,maxTimeStep,isU)
        % select either approximation functions or bisection
        if errs.useApproxFun
            timeStep = priv_approxFun(errs,t,k,k_iter,fullcomp,timeStep,maxTimeStep,isU);

            % approx. function method proposes a value for the time step
            % size which we cannot accept if:
            % 1) value smaller than any lower bound (where errors fulfilled)
            % 2) value larger than any upper bound (where errors not fulfilled)
            % 3) value smaller than any upper bound (where errors fulfilled)
            errs.useApproxFun = ~( ...
                any(timeStep <= [errs.bisect_lb_timeStep_acc, errs.bisect_lb_timeStep_nonacc]) || ...
                any(timeStep >= [errs.bisect_ub_timeStep_acc, errs.bisect_ub_timeStep_nonacc] & ...
                    ~[errs.bisect_ub_accok, errs.bisect_ub_nonaccok]) || ...
                any(timeStep <= [errs.bisect_ub_timeStep_acc, errs.bisect_ub_timeStep_nonacc] & ...
                    [errs.bisect_ub_accok, errs.bisect_ub_nonaccok]) );
        end

        % note: useApproxFun may have changed above!
        if ~errs.useApproxFun
            timeStep = priv_bisection(errs,t,k,fullcomp,timeStep,maxTimeStep,isU);
        end
    end

end


% private methods
methods (Access = private)

    function eps = priv_errOp(S)
        if isa(S,'zonotope')
            % faster computation than calling the interval-constructor
            eps = vecnorm( sum(abs(generators(S)),2) + abs(center(S)) );
        elseif isa(S,'interval')
            eps = vecnorm(max(-infimum(S),supremum(S)));
        else
            % for all others, convert to interval
            I = interval(S);
            eps = radius(or(-I,I));
        end

    end

    function timeStep = priv_approxFun(errs,t,k,k_iter,fullcomp,timeStep,maxTimeStep,isU)
        % predict a time step size which satisfies the error bound
        % e.bound_rem (comprising the remaining error for e.step_acc and
        % e.step_nonacc) using the approximation functions for the
        % individual errors and respecting the maximum admissible time step
        % size due possible changes in the input vector (options.uTransVec)
        
        % timeStep necessary to satisfy eAcc_bound by equating the linearly
        % increasing bound and the approximation function
        
        % in the first step, there is a chance that the initial guess is
        % too large; thus, we use another condition to decrease the initial
        % guess until a reasonable value can be found
        % (note: this works also when e.seq_nonacc{k}(k_iter) = [])
        if k == 1 && any([errs.seq_nonacc{k}(k_iter) / errs.bound_rem{k}(k_iter), ...
                          errs.step_acc{k}(k_iter) / errs.bound_acc{k}(k_iter)] > 1e3)
            timeStep = 0.01 * timeStep;
            return;
        end
        
        % safety factor so that we are more likely to satisfy the bounds
        % in case the approximation functions underestimate the errors
        safetyFactor = 0.90;
        
        % 1. condition: accumulating error needs to satisfy linearly 
        % increasing bound until the time horizon
        % note: if no inputs, then this value becomes Inf
        timeStep_pred_linBound = safetyFactor * ...
            (1/errs.coeff_PUtkplus1_a * ...
            (errs.bound_remacc(k) / (errs.tFinal - t) - errs.coeff_PUtkplus1_b));
        
        timeStep_pred_erem = Inf;
        if fullcomp
            % 2. condition: both errors need to satisfy erem for the
            % current step ...this yields a quadratic equation
            temp_a = errs.coeff_PUtkplus1_a + errs.coeff_F_a + errs.coeff_G_a;
            temp_b = errs.coeff_PUtkplus1_b + errs.coeff_F_b + errs.coeff_G_b ...
                        + errs.coeff_linComb_b + errs.coeff_PUtauk_b;
            temp_c = -errs.bound_rem{k}(k_iter);
            
            % compute two solutions of quadratic equation (one is < 0)
            % choose maximum of predicted timeSteps -> positive solution
            timeStep_pred_erem = safetyFactor * max(...
                [1/(2*temp_a) * (-temp_b + sqrt(temp_b^2 - 4*temp_a*temp_c)); ...
                 1/(2*temp_a) * (-temp_b - sqrt(temp_b^2 - 4*temp_a*temp_c))] );
        end
        
        % update timeStep by minimum of predicted timeSteps so that both
        % error bounds are likely to be satisfied; also, the maximum
        % admissible for the current step needs to be respected, as the
        % input vector has to be constant over any single step
        timeStep = min([maxTimeStep;timeStep_pred_linBound;timeStep_pred_erem]);
    end

    function timeStep = priv_bisection(errs,t,k,fullcomp,timeStep,maxTimeStep,isU)

        % special handling for first step
        if k == 1
            eacctotal = 0;
        else
            eacctotal = errs.cum_acc(k-1);
        end
        
        % 1. acc errors
        timeStep_pred_eacc = Inf;
        slope_acc = 0;
        if isU
        if errs.bisect_ub_accok
            % extrapolate
            errorbound_accend_tFinal = (errs.emax - errs.bound_red_max) / errs.tFinal;
            slope_acc = (errs.bisect_ub_acc - errs.bisect_lb_acc) / ...
                (errs.bisect_ub_timeStep_acc - errs.bisect_lb_timeStep_acc);
            
            if slope_acc < errorbound_accend_tFinal
                % slope of error bound larger than slope of error (this occurs due
                % to assumed linearity... actually, the slope of the error will
                % increase for larger values due to higher-order terms)
                timeStep_pred_eacc = Inf;
            else
                timeStep_add = ( errorbound_accend_tFinal * (t + errs.bisect_lb_timeStep_acc) ...
                    - eacctotal - errs.bisect_lb_acc ) / (slope_acc - errorbound_accend_tFinal);
                timeStep_pred_eacc = errs.bisect_lb_timeStep_acc + timeStep_add;
            end
        else
            % bisection
            if k == 1 % constant in first step as differences very large
                factor = 0.5;
            else % weighted, using information that decrease at least quadratically
                % ...estimate so that bound fulfilled at 95%
                factor = sqrt( (0.95 - errs.bisect_lb_acc_perc) / ...
                    (errs.bisect_ub_acc_perc - errs.bisect_lb_acc_perc) );
            end
            timeStep_pred_eacc = (errs.bisect_lb_timeStep_acc + ...
                factor * (errs.bisect_ub_timeStep_acc - errs.bisect_lb_timeStep_acc));
        end
        end
        
        % 2. nonacc errors
        timeStep_pred_enonacc = Inf;
        if fullcomp
        if errs.bisect_ub_nonaccok
            % extrapolate
            slope_nonacc = (errs.bisect_ub_nonacc - errs.bisect_lb_nonacc) / ...
    	        (errs.bisect_ub_timeStep_nonacc - errs.bisect_lb_timeStep_nonacc);
            if isU
                timeStep_add = (errs.emax - errs.bound_red_max*t/errs.tFinal - eacctotal ...
                    - errs.bisect_lb_acc - errs.bisect_lb_nonacc ) / ...
                    ( slope_acc + slope_nonacc - errs.bound_red_max/errs.tFinal );
            else 
                % shorter way
                timeStep_add = (errs.emax - errs.bisect_lb_nonacc) / slope_nonacc;
            end
            timeStep_pred_enonacc = errs.bisect_lb_timeStep_nonacc + timeStep_add;
        else
            % bisection
            if k == 1 % constant in first step as differences very large
                factor = 0.5;
            else % weighted, estimate so that bound fulfilled at 95%
                factor = (0.95 - errs.bisect_lb_nonacc_perc) / ...
                    (errs.bisect_ub_nonacc_perc - errs.bisect_lb_nonacc_perc);
            end
            timeStep_pred_enonacc = (errs.bisect_lb_timeStep_nonacc + ...
                    factor * (errs.bisect_ub_timeStep_nonacc - errs.bisect_lb_timeStep_nonacc));
        end
        % for safety... (if slopes misbehave)
        if timeStep_pred_enonacc < 0
            timeStep_pred_enonacc = Inf;
        end
        end
        
        % find minimum
        timeStep = min([timeStep_pred_enonacc,timeStep_pred_eacc,maxTimeStep]);
    end
end


% print/plot functions
methods
    
    function print(errs)
        % print information on command window: to be called after all
        % computations have been done
        
        % init table
        hvalues = {'Step','Time Step Size',...
            'Non-acc. Error','Non-acc. Error Bound',...
            'Acc. Error','Acc. Error Bound',...
            'Red. Error','Red. Error Bound'};
        formats = {'i','.4e','.3e','.3e','.3e','.3e','.3e','.3e'};
        table = CORAtable('single',hvalues,formats);
        
        % print table
        table.printHeader()
        for k=1:numel(errs.timeSteps)
            table.printContentRow({...
                k,errs.timeSteps{k},...
                errs.seq_nonacc{k},errs.bound_rem{k},...
                errs.step_acc{k},errs.bound_acc{k},...
                errs.step_red(k),errs.bound_red{k},...
                });
        end
        table.printFooter();
        
    end

    function plot(errs)
        % plot non-accumulating, accumulating, reduction, and total error,
        % including the corresponding bounds

        % staircase function for time
        t = cumsum(cell2mat(errs.timeSteps));
        t_plot = [0;repelem(t(1:end-1),2);t(end)];

        figure;
        % non-accumulating errors and bound
        subplot(2,2,1); hold on; box on;
        e_plot = repelem(cell2mat(errs.seq_nonacc),2);
        h_err = plot(t_plot,e_plot);
        e_bound_plot = repelem(cell2mat(errs.bound_rem),2);
        h_bound = plot(t_plot,e_bound_plot);
        xlabel("Time"); ylabel("Non-accumulating error (bound)");
        axis([0,errs.tFinal,0,1.1*max([e_plot;e_bound_plot])]);
        legend([h_err,h_bound],'Error','Error bound');

        % accumulating errors and bound
        subplot(2,2,2); hold on; box on;
        e_plot = repelem(cell2mat(errs.step_acc),2);
        h_err = plot(t_plot,e_plot);
        e_bound_plot = repelem(cell2mat(errs.bound_acc),2);
        h_bound = plot(t_plot,e_bound_plot);
        xlabel("Time"); ylabel("Accumulating error (bound)");
        axis([0,errs.tFinal,0,1.1*max([e_plot;e_bound_plot])]);
        legend([h_err,h_bound],'Error','Error bound');

        % reduction error and bound
        subplot(2,2,3); hold on; box on;
        e_plot = repelem(errs.step_red,2);
        h_err = plot(t_plot,e_plot);
        e_bound_plot = repelem(cell2mat(errs.bound_red),2);
        h_bound = plot(t_plot,e_bound_plot);
        xlabel("Time"); ylabel("Reduction error (bound)");
        axis([0,errs.tFinal,0,1.1*max([e_plot;e_bound_plot])]);
        legend([h_err,h_bound],'Error','Error bound');

        % total error
        subplot(2,2,4); hold on; box on;
        e_plot = repelem(cell2mat(errs.seq_nonacc),2) + ...
            repelem(cell2mat(errs.step_acc),2) + repelem(errs.step_red,2);
        h_err = plot(t_plot,e_plot);
        h_bound = plot([0,errs.tFinal],[errs.emax,errs.emax]);
        xlabel("Time"); ylabel("Total error (bound)");
        axis([0,errs.tFinal,0,1.1*max([e_plot;errs.emax])]);
        legend([h_err,h_bound],'Error','Error bound');
        
    end

end

end

% ------------------------------ END OF CODE ------------------------------
