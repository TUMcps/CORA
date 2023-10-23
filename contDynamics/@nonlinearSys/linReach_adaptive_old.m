function [Rti,Rtp,options] = linReach_adaptive_old(obj,options,Rstart)
% linReach_adaptive_old - computes the reachable set after linearization
%
% Syntax:
%    [Rti,Rtp,options] = linReach_adaptive_old(obj,options,Rstart)
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
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Mark Wetzlinger
% Written:       11-December-2020
% Last update:   15-January-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set linearization error
error_adm = zeros(obj.dim,1);

% time step size + set abstraction order ----------------------------------
% finite horizon depending on varphi, init at start, control below
lastStep = false; veryfirststep = false; timeStepequalDelta = false;
if options.i == 1
    veryfirststep = true;
    options.timeStep = options.tFinal * 0.01; % initial guess
    finitehorizon = options.timeStep;
elseif options.timeStep > options.tFinal - options.t
    lastStep = true;
    options.timeStep = options.tFinal - options.t;
    finitehorizon = options.timeStep;
elseif options.i > 1
    % take last one and estimate so that varphi = zetaphi
    finitehorizon = options.finitehorizon(options.i-1) * ...
        (options.decrFactor - options.zetaphi) / ...
        (options.decrFactor - options.varphi(options.i-1));
    % cap for change in finite horizon
    options.timeStep = finitehorizon;
    if finitehorizon > options.tFinal - options.t
        options.timeStep = options.tFinal - options.t;
        finitehorizon = options.timeStep;
    end
end

% iteration counter run:
% ...0 only for first step in lin
% ...otherwise: 1: compute finitehorizon (use last deltat for varphi), 2: compute optimal deltat
options.run = 0;
% adaptive tensor order (first step: pre-computation by comparison on start sets)
if veryfirststep && strcmp(options.alg,'lin')
    % compute L for start set to decide abstraction order
    options = aux_lbL(obj,options,Rstart);
end
thistO = options.tensorOrder;
options.run = options.run + 1;
% -------------------------------------------------------------------------
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
    
    perfIndCounter = 1; perfInds = []; Lconverged = true;
    while true
        abscount = abscount + 1;
        % estimate the abstraction error
        appliedError = 1.1*error_adm; % error_adm always interval
        if any(error_adm)
            % set resulting from abstraction error
            errG = diag(appliedError);
            Verror = zonotope([0*appliedError,errG(:,any(errG,1))]);
            % compute abstraction error as input solution
            RallError = errorSolution_adaptive(linSys,options,Verror);
        else
            % only used in options.i == 1, but important for adaptive
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
                throw(CORAerror('CORA:wrongFieldValue','options.alg',...
                    {'lin','poly'}));
            end
        catch ME
            if strcmp(ME.identifier,'reach:setoutofdomain')
                Lconverged = false; break
            else
                rethrow(ME);
            end
        end
        
        % compare linearization error with the maximum admissible error
        perfIndCurr = max(trueError ./ appliedError);
        if perfIndCurr <= 1 || ~any(trueError)
            % error converged: compute tensorOrder for next step (only lin)
            if strcmp(options.alg,'lin') && options.run == 2
                % compute L for other tensorOrder (2/3)
                if options.tensorOrder == 2
                    L0_2 = trueError; options.tensorOrder = 3;
                    [~,~,L0_3,~] = abstractionError_adaptive(obj,options,Rmax,Rtemp);
                else
                    L0_3 = trueError; options.tensorOrder = 2;
                    [~,~,L0_2,~] = abstractionError_adaptive(obj,options,Rmax,Rtemp);
                end

                % cut linear dimensions
                L0_2lin = L0_2(L0_2 ~= 0); L0_3lin = L0_3(L0_3 ~= 0);

                % compare Lagrange remainder to decide which
                %    tensorOrder to use (for next step)
                if all( L0_3lin ./ L0_2lin > options.zetaK ); nexttO = 2;
                else;                                         nexttO = 3;
                end
                options.tensorOrder = thistO; % for safety
            end
            break
        elseif perfIndCounter > 1
            perfInds(perfIndCounter) = perfIndCurr;
            if perfIndCounter > 2 && perfInds(perfIndCounter) > perfInds(perfIndCounter-1)
                Lconverged = false; break
            end
        end
        
        error_adm = trueError;
        perfIndCounter = perfIndCounter + 1;
    end
    
    % if either L was out of domain or did not converge
    if ~Lconverged
        options.timeStep = options.timeStep * 0.5;
        finitehorizon = options.timeStep;
        error_adm = zeros(obj.dim,1);
        continue
    end
    % ... now containment of L ensured
    
    if veryfirststep
        % very first run through here
        options = aux_getPowers(linSys,options,VerrorDyn);
    end
    
    % compute set of abstraction errors
    [Rerror,options] = errorSolution_adaptive(linSys,options,VerrorDyn,VerrorStat);
    
    % measure abstraction error
    if isa(Rerror,'polyZonotope')
        abstrerr = sum(abs(Rerror.GI),2)';
    else
        abstrerr = sum(abs(generators(Rerror)),2)';
    end
    
    % two runs in each step
    if options.run == 1 && ~lastStep
        % ... run using finite horizon
        if options.i == 1
            % first step more cumbersome as there is no previous knowledge
            if veryfirststep
                % shrink time step by decrFactor to get varphi (only options.i == 1)
                veryfirststep = false;
                abstrerr_Delta = abstrerr;
                Rerror_Delta = Rerror;
                Rti_Delta = Rlinti; Rtp_Delta = Rlintp;
                linx_Delta = obj.linError.p.x;
                zetaP_Delta = zetaP;
                options.timeStep = options.decrFactor * options.timeStep;
                error_adm = zeros(obj.dim,1);
                continue;
            end
            % take worst-case current gain (use same dimensions as for lin)
            varphi = max( abstrerr(options.lindims) ./ abstrerr_Delta(options.lindims) );
            if options.i == 1 && varphi < options.zetaphi
                % varphi not large enough ... decrease again
                finitehorizon = options.timeStep;
                abstrerr_Delta = abstrerr;
                Rerror_Delta = Rerror;
                Rti_Delta = Rlinti; Rtp_Delta = Rlintp;
                linx_Delta = obj.linError.p.x;
                zetaP_Delta = zetaP;
                options.timeStep = options.decrFactor * options.timeStep;
                error_adm = zeros(obj.dim,1);
                continue;
            end
            options.finitehorizon(options.i,1) = finitehorizon;
            % estimate near-optimal delta t from optimization function
            [options.timeStep,~] = aux_optimaldeltat(Rstart,Rerror_Delta,...
                finitehorizon,varphi,zetaP_Delta,options);
            % reuse information (scrap once pre-computation replaced)
            error_adm = zeros(obj.dim,1);
            
        elseif options.i > 1
            abstrerr_Delta = abstrerr;
            Rerror_Delta = Rerror;
            Rti_Delta = Rlinti; Rtp_Delta = Rlintp;
            linx_Delta = obj.linError.p.x;
            [options.timeStep,~] = aux_optimaldeltat(Rstart,Rerror_Delta,...
                finitehorizon,options.varphi(options.i-1),zetaP,options);
            error_adm = zeros(obj.dim,1);
        end
        
    elseif options.run == 2 || lastStep
        % ... run using tuned time step size / scaled time step size for varphi
        if timeStepequalDelta
            options.varphi(options.i,1) = ...
                max( abstrerr(options.lindims) ./ abstrerr_Delta(options.lindims) );
        elseif ~lastStep
            options.varphi(options.i,1) = aux_varphiest(finitehorizon,options.timeStep,...
                Rerror_Delta,Rerror,options.decrFactor);
        end
        options.finitehorizon(options.i,1) = finitehorizon;
        % second run finished -> exit
        options.abscount(options.i,1) = abscount;
        break
    end
    
    % update counter
    options.run = options.run + 1;
    
    if options.timeStep == finitehorizon
        % required Rerror and Rlin already there, but compute
        % decrFactor*finitehorizon for usable varphi
        timeStepequalDelta = true;
        options.timeStep = options.timeStep * options.decrFactor;
    end 
    
end

% use correct Rerror in case optimal time step size is equal to finitehorizon
if timeStepequalDelta
    % reset time step size for time increment in loop outside
    options.timeStep = finitehorizon;
    % choose sets computed by Delta (last run: finitehorizon*decrFactor)
    Rti = Rti_Delta + linx_Delta; Rtp = Rtp_Delta + linx_Delta;
    Rerror = Rerror_Delta;
else
    % translate reachable sets by linearization point
    Rti = Rlinti + obj.linError.p.x; Rtp = Rlintp + obj.linError.p.x;
end

% add the abstraction error to the reachable sets
if isa(Rerror,'polyZonotope')
    Rti=exactPlus(Rti,Rerror); Rtp=exactPlus(Rtp,Rerror);
else
    Rti=Rti+Rerror; Rtp=Rtp+Rerror;
end

% varphi prediction / regression of next step (change once all ready)
if strcmp(options.alg,'lin') && ~lastStep && nexttO ~= thistO
    % switch of tensorOrder: compute also other Rerror
    % (required for varphi of next step)
    options.tensorOrder = nexttO;
end

end


% Auxiliary functions -----------------------------------------------------

function varphi = aux_varphiest(Delta,deltatlast,Rerr,Rerrlast,decrFactor)

% compute set sizes
if isa(Rerr,'zonotope') && isa(Rerrlast,'zonotope')
    rerr1 = vecnorm(sum(abs(generators(Rerr)),2),2);
    rerrk = vecnorm(sum(abs(generators(Rerrlast)),2),2);
else
    rerr1 = vecnorm(sum(abs(Rerr.GI),2),2);
    rerrk = vecnorm(sum(abs(Rerrlast.GI),2),2);
end

k = Delta / deltatlast;
rhs = k * rerrk / rerr1;
kprime = - log(k) / log(decrFactor);
kprime_floor = floor(kprime);
kprime_ceil = ceil(kprime);

tol = 0.001;

% starting lb / ub by 0 and decrFactor
varphi_floor_lb = 0;
varphi_ceil_lb = 0;
varphi_floor_ub = decrFactor;
varphi_ceil_ub = decrFactor;

varphi_floor = varphi_floor_lb;
varphi_ceil = varphi_ceil_lb;

% binary search
counter = 0;
while true
    allvarphi_floor = ( varphi_floor + (1 - decrFactor .^ [0:kprime_floor]) ...
        * (decrFactor - varphi_floor) ) / decrFactor;
    allvarphi_ceil = ( varphi_ceil + (1 - decrFactor .^ [0:kprime_ceil]) ...
        * (decrFactor - varphi_ceil) ) / decrFactor;
    res_floor = prod(allvarphi_floor) - rhs;
    res_ceil = prod(allvarphi_ceil) - rhs;
    
    % termination condition
    if abs(res_floor) < tol && abs(res_ceil) < tol
        break
    end
    if res_floor < 0
        varphi_floor_lb = varphi_floor;
        varphi_floor = 0.5 * (varphi_floor + varphi_floor_ub);
    else
        varphi_floor_ub = varphi_floor;
        varphi_floor = 0.5 * (varphi_floor_lb + varphi_floor);
    end
    if res_ceil < 0
        varphi_ceil_lb = varphi_ceil;
        varphi_ceil = 0.5 * (varphi_ceil + varphi_ceil_ub);
    else
        varphi_ceil_ub = varphi_ceil;
        varphi_ceil = 0.5 * (varphi_ceil_lb + varphi_ceil);
    end
    counter = counter + 1;
    if counter > 10000
        throw(CORAerror('CORA:notConverged','Adaptive tuning'));
    end
end
varphi = varphi_floor + (kprime - floor(kprime)) * (varphi_ceil - varphi_floor);


end

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

rabskest = rerr1 / k(bestIdxnew) * varphiprod(bestIdxnew);

end

function options = aux_lbL(obj,options,Rstartset)
% computes the lower bound for the Lagrange remainder
% ... that is, L and Rerr for the start set of the current step

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

function options = aux_getPowers(obj,options,Vdyn)
% computes linear exponents for convergence of abstraction error
% only called once!

n = obj.dim;
Apower = eye(n);
Vgens = generators(Vdyn);
powers = ones(n,1);

% powers with linear gain
ApowerV = Apower * Vgens;
powerIncr = ~any(ApowerV,2);
powers = powers + ones(n,1) .* powerIncr;

options.lindims = powers == 1;

end

% ------------------------------ END OF CODE ------------------------------
