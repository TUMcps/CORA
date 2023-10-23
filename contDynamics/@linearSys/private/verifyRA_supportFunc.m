function [res,fals,savedata] = verifyRA_supportFunc(obj,params,options,spec)
% verifyRA_supportFunc - verification for linear systems
%    via reach-avoid with support functions:
%    quicker verification algorithm based on non-trivial
%    exploitations of the standard propagation-based wrapping-free
%    reachability algorithm to the extent that reachable sets only computed
%    implicitly with respect to their distance to an unsafe set
%    caution: specification needs to be given as a safe set!
%
% Syntax:
%    res = verifyRA_supportFunc(obj,params,options,spec)
%
% Inputs:
%    obj - linearSys object
%    params - model parameters
%    options - algorithm parameters
%    spec - safe set as object of class specification
%
% Outputs:
%    res - boolean (true if specification verified, false if not)
%    fals - struct containing falsifying trajectory
%           .x0 ... point from initial set
%           .u  ... piecewise-constant input values
%           .tu ... switching times of .u
%    savedata - distances for plotting
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/verify

% Authors:       Mark Wetzlinger
% Written:       04-April-2022
% Last update:   22-April-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% integrate in verify?
options = validateOptions(obj,mfilename,params,options);
% remove after development phase
options.verbose = false;

% start stopwatch
onlyLoop = tic;
% initializations ---------------------------------------------------------
if options.verbose
    disp("...initializations");
end

% init time shift
options.tFinal = options.tFinal - options.tStart;

% initialize satisfaction of specification
res = false;

% get safe sets and unsafe sets G from specification
% TODO: remove once integration in verify has been done
[safeSet,unsafeSet] = aux_getSetsFromSpec(spec);
% get halfspaces from specifications (transform all to safe sets)
nrSpecs = length(safeSet) + length(unsafeSet);
Cs = zeros(nrSpecs,obj.nrOfOutputs);
ds = zeros(nrSpecs,1);
for i=1:length(safeSet)
    % if all values satisfy C*x <= d, system is safe
    Cs(i,:) = safeSet{i}.set.A;
    ds(i) = safeSet{i}.set.b;
end
for j=1:length(unsafeSet)
    % if any value satisfies C*x <= d, system is unsafe
    % -> rewrite to safe set: C*x >= d  <=>  -C*x <= -d
    Cs(i+j,:) = -unsafeSet{j}.set.A;
    ds(i+j) = -unsafeSet{j}.set.b;
end

% time intervals where specifications are not yet verified
% TODO: integrate in validateOptions later...
specUnsat = cell(nrSpecs,1);
for s=1:nrSpecs
    specUnsat{s} = [0,options.tFinal];
end

% initialize all variables for exponential (propagation) matrix
expmat = aux_initExpmat(obj);

% rewrite equations in canonical form
[options,origInput] = aux_canonicalForm(obj,options);

% compute directions for support function evaluation: the specifications
%    Cs*y <= ds
% are defined in the output space, thus we insert y=C*x to obtain
%    Cs*y <= d  ->  Cs*(C*x) <= d  ->  (Cs * C)*x <= d
% and map all Cs to l(i) = Cs(i)*C for further use
l = (Cs * full(obj.C))';

% scaling factor due to output matrix (only used for heuristic to adapt the
% time step size)
options.normC = norm(full(obj.C));

% init all auxiliary structs and flags
[dist,fals,isU,G_U,isu,isuconst,options] = aux_initStructsFlags(obj,options);

% initialize time step size as largest whole number divisor over all
% individual piecewise-constant inputs
nrStepsStart = 100;
timeStep = aux_timeStep(options.tFinal/nrStepsStart,options.tFinal,isuconst,options.tu);
% maximum factor by which time step size may decrease in one iteration
timeStepFactor = 1/50;
% factor by which time step size decreases if a sum does not converge
timeStepFactor_nonconverged = 0.2;
% method from paper: fixed factor
timeStepFactor_fixed = 0.2;

% init savedata
savedata.tFinal = options.tFinal;
savedata.isU = isU;
savedata.nrSpecs = nrSpecs;
savedata.iterations = 0;

% plots on (only debugging) or off
isplot = false; blobsize = 3;
% -------------------------------------------------------------------------

% quick check: does start set already intersect the unsafe set?

% debug ---
if isplot
    close all;
    figure; 
    for s=1:nrSpecs
        subplot(1,nrSpecs,s); hold on; box on;
        plot([0,options.tFinal],[ds(s) ds(s)],'Color',colorblind('r'));
        scatter(0,supremum(interval(l(:,s)' * options.R0 + options.V ...
            + options.vTransVec(:,1))),blobsize,'k','filled');
    end
end
% debug ---

% shift distance to halfspace by influence from additive output set V
% (from y = C*x + v, v \in V... only accounting for constant V!),
% so that now we have already dealt with V for all sets;
% note: we use only Cs instead of l = Cs*C as V is defined in output space!
ds = ds - Cs*center(options.V) - sum(abs(Cs*generators(options.V)),2);

% read data from initial set (these remain constant for the algorithm!)
c_X0 = center(options.R0);
G_X0 = generators(options.R0);

% compute distance of initial output set to each specification
dist_affine_tp_0 = (-ds + l'*c_X0 + sum(abs(l'*G_X0),2) + Cs*options.vTransVec(:,1))';

% check initial output set (outside of verification loop below)
if any(dist_affine_tp_0 > 0)
    % check which spec causes falsification
    falsIdx = find(dist_affine_tp_0 > 0,1,'first');
    % time where falsification occurs: initial time
    fals.tFinal = 0;
    % most critical initial state (via support function evaluation)
    fals.x0 = c_X0 + G_X0*sign(l(:,falsIdx)'*G_X0)';
    % no input trajectory needed: same as input trajectory = 0 at all times
    fals.u = zeros(dim(origInput.U),1);
    fals.tu = 0;
    % measure elapsed time
    savedata.tComp = toc(onlyLoop);
    % provably falsified -> exit!
    return;
end


% iteration counter
savedata.iterations = 0;
% main loop: continues until specification is verified/falsified
while true
    
    % increment counter
    savedata.iterations = savedata.iterations + 1;

    % compute number of steps (we use fixed time step sizes)
    nrSteps = round(options.tFinal/timeStep);

    % save data
    savedata.timeStep = timeStep;
    savedata.nrSteps = nrSteps;

    % log
    if options.verbose
        disp("Iteration " + savedata.iterations + ": no. of steps = " + nrSteps + ...
            " (time step size = " + timeStep + ", horizon = " + options.tFinal + ")");
    end

    % initialize contributions to distance to unsafe set from piecewise
    % constant input solution, affine solution, full curvature
    dist.Pu = zeros(nrSteps+1,nrSpecs);
    dist.affine_tp = [dist_affine_tp_0; NaN(nrSteps,nrSpecs)];
    dist.Cbloat = zeros(nrSteps,nrSpecs);
    % initialize contribution to distance to unsafe set from input vector
    % in output equation (assuming constant input vector for now,
    % recomputed below if necessary); note: included in dist.affine_tp!
    dist.vTrans = repmat((Cs*options.vTransVec(:,1))',nrSteps+1,1);

    % compute exponential matrix for given time step size (+ transposed)
    if options.verbose
        disp("...compute propagation matrix");
    end
    expmat.Deltatk = expm(obj.A * timeStep);
    expmat.Deltatk_T = expmat.Deltatk';

    % compute constant input solution
    if isu
        % assign constant shifts
        [u,vnext,tnextSwitch] = aux_uTrans_vTrans(0,timeStep,...
            options.uTransVec,options.vTransVec,options.tu,0,...
            options.uTransVec(:,1),options.vTransVec(:,1));
        Pu = aux_Pu(obj,u,expmat,timeStep);
        % initialize accumulated solution
        Pu_total = zeros(obj.dim,1);
        % reduce time step size if computation has not converged
        if ~expmat.conv
            timeStep = timeStep * timeStepFactor_nonconverged;
            expmat.conv = true;
            if options.verbose
                disp("...reduce time step size");
            end
            continue;
        end
    end

    % compute interval matrices required for curvature error sets
    if options.verbose
        disp("...compute interval matrices for curvature errors");
    end
    [expmat,intmatF,intmatG] = aux_intmat(obj,isu,expmat,timeStep);
    % reduce time step size if computation has not converged
    if ~expmat.conv
        timeStep = timeStep * timeStepFactor_nonconverged;
        expmat.conv = true;
        if options.verbose
            disp("...reduce time step size");
        end
        continue;
    end

    % original computation (not taking time-varying u into account)
%     Cbloat = intmatF.set * options.R0 + intmatG.set * u;
%     c_Cbloat = center(Cbloat); G_Cbloat = generators(Cbloat);
    % speed-up
    cx_Cbloat = intmatF.center * c_X0;
    G1x_Cbloat = intmatF.center * G_X0;
    G2x_Cbloat = intmatF.rad * sum(abs([c_X0 G_X0]),2);
    if isu
        cu_Cbloat = intmatG.center * u;
        Gu_Cbloat = intmatG.rad * abs(u);
    end
    % note: G2x and Gu require diag(.), but we can omit it here in order to
    % speed up the computation of dist.Cbloat!

    % first curvature distance out of the loop (better indexing)
    dist.Cbloat(1,:) = l'*cx_Cbloat + sum(abs(l'*G1x_Cbloat),2) ...
        + abs(l'*G2x_Cbloat);
    if isu
        dist.Cbloat(1,:) = dist.Cbloat(1,:)' ...
            + l'*cu_Cbloat + abs(l'*Gu_Cbloat);
    end

    % initialize variables for back-propagated direction of support
    % function (including transpose)
    l_prop = cell(nrSpecs,1);
    l_prop_T = cell(nrSpecs,1);
    for s=1:nrSpecs
        l_prop{s} = [l(:,s), zeros(obj.dim,nrSteps)];
        l_prop_T{s} = l_prop{s}';
    end

    % log
    if options.verbose
        disp("...compute affine solutions"); 
    end
    if nrSteps > 10
        logVar = 0.1;
    end

    % compute distance of affine time-point solutions H(tk) to unsafe set
    for k=1:nrSteps
        % log
        if options.verbose && nrSteps > 10 && k > round(nrSteps * logVar)
            if logVar == 0.1
                fprintf('...');
            end
            fprintf([num2str(100*logVar), '%%, ']);
            logVar = logVar + 0.1;
        end

        % back-propagate direction of support function (incl. transpose)
        for s=1:nrSpecs
            l_prop{s}(:,k+1) = expmat.Deltatk_T * l_prop{s}(:,k);
            l_prop_T{s}(k+1,:) = l_prop{s}(:,k+1)';
        end

        % propagate constant input solution
        if isu
            if ~isuconst
                % compute next(!) Pu and contribution to bloating term
                [u,vnext,tnextSwitch] = aux_uTrans_vTrans((k-1)*timeStep,...
                    timeStep,options.uTransVec,options.vTransVec,...
                    options.tu,tnextSwitch,u,vnext);

                if k > 1 && ~all(withinTol(u,u_prev))
                    % only update Pu if necessary
                    Pu = aux_Pu(obj,u,expmat,timeStep);
                    % update values for Cbloat
                    cu_Cbloat = intmatG.center * u;
                    % note: Gu requires diag(.), but we can omit it here in
                    % order to speed up the computation of dist.Cbloat!
                    Gu_Cbloat = intmatG.rad * abs(u);
                end

                % store u for comparison to avoid re-computations
                u_prev = u;
            end
            
            % accumulate particular solution
            Pu_total = expmat.Deltatk * Pu_total + Pu;

            % actual distance computation (accumulates over time as the
            % particular constant input solution does too)
            dist.Pu(k+1,:) = (l' * Pu_total)';

            % distance contribution from input vector in output equation
            dist.vTrans(k+1,:) = (Cs*vnext)';
        end

        % propagate time-point solution H(tk), compute distance of output
        % time-point solution Y(tk) to unsafe set
        for s=1:nrSpecs
            dist.affine_tp(k+1,s) = -ds(s) ...
                + l_prop_T{s}(k+1,:) * c_X0 ...
                + sum(abs(l_prop_T{s}(k+1,:)*G_X0),2) ...
                + dist.Pu(k+1,s) ...
                + dist.vTrans(k+1,s);
        end

        % in case of an intersection, already falsified
        if any(dist.affine_tp(k+1) > 0)
            % provably falsified!
            savedata.tComp = toc(onlyLoop);

            % log
            if options.verbose
                fprintf('...falsification detected!\n');
            end

            % debug ---
            if isplot
                for s=1:nrSpecs
                    subplot(1,nrSpecs,s); hold on;
                    scatter(timeStep*(0:nrSteps),...
                        ds(s)+dist.affine_tp(:,s),blobsize,'k','filled');
                end
            end
            % debug ---

            % check which spec causes falsification
            falsIdx = find(dist.affine_tp(k+1) > 0,1,'first');
            % time where falsification occurs
            fals.tFinal = timeStep * k;
            % starting point x0 which yields the falsifying trajectory
            fals.x0 = c_X0 + G_X0*sign(l_prop_T{falsIdx}(k+1,:)*G_X0)';
            % falsifying piecewise-constant input is equal to uTransVec
            % until current point in time...
            if ~isu || isuconst
                fals.u = origInput.u;
                fals.tu = 0;
            else
                idx = find(k*timeStep >= options.tu,1,'last');
                fals.u = origInput.u(:,1:idx);
                fals.tu = options.tu(1:idx);
            end
            % save data
            savedata.fals_tFinal = fals.tFinal;
            savedata.dist_affinetp = dist.affine_tp;
            return;
        end

        % computation of distance contribution from curvature term and
        % offset in output equation required for the time-interval affine
        % solution H(tauk)
        if k < nrSteps
            for s=1:nrSpecs
                dist.Cbloat(k+1,s) = l_prop_T{s}(k+1,:)*cx_Cbloat ...
                    + sum(abs(l_prop_T{s}(k+1,:)*G1x_Cbloat)) ...
                    + abs(l_prop_T{s}(k+1,:)*G2x_Cbloat);
                if isu
                    dist.Cbloat(k+1,s) = dist.Cbloat(k+1,s) ...
                        + l_prop_T{s}(k+1,:)*cu_Cbloat ...
                        + abs(l_prop_T{s}(k+1,:)*Gu_Cbloat);
                end
            end
        end

    end
    % end log
    if options.verbose && nrSteps > 10
        fprintf('100%%\n');
    end

    % debug ---
    if isplot
        for s=1:nrSpecs
            subplot(1,nrSpecs,s); hold on;
            scatter(timeStep*(0:nrSteps),...
                ds(s)+dist.affine_tp(:,s),blobsize,'k','filled');
        end
    end
    % debug ---
    
    % compute distance for time-interval solution H(tauk): add curvature
    % to both start and end set of each H(tk) -> 2 columns (note: the
    % piecewise-constant input vector vTransVec is included in
    % dist.affine_tp)
    dist.affine_ti = cell(nrSpecs,1);
    for s=1:nrSpecs
        dist.affine_ti{s} = [dist.affine_tp(1:end-1,s),...
            dist.affine_tp(2:end,s)] + dist.Cbloat(:,s);
    end
    % TODO: inner-approx of affine_ti with -dist.Cbloat!
    
    % debug ---
    if isplot
        for s=1:nrSpecs
            subplot(1,nrSpecs,s); hold on;
            plot(timeStep * [0;repelem(1:nrSteps-1,1,2)';nrSteps],...
                reshape((ds(s)+dist.affine_ti{s})',nrSteps*2,1),...
                'Color',colorblind('b'));
        end
    end
    % debug ---

    % gather already verified time intervals (unless input set provided)
    if ~isU
        if options.verbose
            disp("...check distances");
        end

        % quick version
        if all(cellfun(@(x) all(all(x<0)),dist.affine_ti,'UniformOutput',true))
            specUnsat = cell(nrSpecs,1);
        elseif any(cellfun(@(x) any(any(x<0)),dist.affine_ti,'UniformOutput',true))
            for s=1:nrSpecs
                t = [0,timeStep];
                for k=1:nrSteps
                    if all(dist.affine_ti{s}(k,:) < 0)
                        specUnsat{s} = aux_removeFromUnsat(specUnsat{s},t);
                    end
                    t = t + timeStep;
                end
            end
        end

        % suggest refinement based on H(tauk) (important if no input set)
        if ~all(cellfun(@isempty,specUnsat,'UniformOutput',true))
            % shorten time horizon if possible
            options.tFinal = max(cellfun(@(x) x(end,2),specUnsat,'UniformOutput',true));
            % method from paper: fixed factor
            timeStep = timeStep * timeStepFactor_fixed;
            % propose smaller time step size
%             timeStep = aux_affineTimeStep(obj,dist,options.tFinal,...
%                 timeStep,timeStepFactor,isuconst,options.tu,options.normC);
        end


    % remainder only required if there is an input to the system
    else
        if options.verbose
            disp("...compute inner/outer-approximation of particular solution");
        end
        % compute inner-approximation underPU(Delta t)
        underPU_G = aux_underPU(obj,G_U,expmat,timeStep);
        % reduce time step size if computation has not converged
        if ~expmat.conv
            timeStep = timeStep * timeStepFactor_nonconverged;
            expmat.conv = true;
            if options.verbose
                disp("...reduce time step size");
            end
            continue;
        end
        % compute outer-approximation overPU(Delta t)
        overPU_G = aux_overPU(obj,G_U,expmat,timeStep);
        % reduce time step size if computation has not converged
        if ~expmat.conv
            timeStep = timeStep * timeStepFactor_nonconverged;
            expmat.conv = true;
            if options.verbose
                disp("...reduce time step size");
            end
            continue;
        end

        % propagate distance shifts (index k for solution at tk since both
        % inner- and outer-approximation are 0 at time t0)
        dist.underPU = zeros(nrSteps+1,nrSpecs);
        dist.overPU = zeros(nrSteps+1,nrSpecs);
        for s=1:nrSpecs
            dist.underPU(:,s) = [0;cumsum(sum(abs(l_prop_T{s}(1:end-1,:)*underPU_G),2))];
            dist.overPU(:,s) = [0;cumsum(sum(abs(l_prop_T{s}(1:end-1,:)*overPU_G),2))];
        end

        % debug ---
        if isplot
            for s=1:nrSpecs
                subplot(1,nrSpecs,s); hold on;
                % inner-approximation
                scatter(timeStep*(0:nrSteps),...
                    ds(s)+dist.affine_tp(:,s)+dist.underPU(:,s),...
                    blobsize,'g','filled');
                % outer-approximation
                plot(timeStep * [0;repelem(1:nrSteps-1,1,2)';nrSteps],...
                    reshape((ds(s)+dist.affine_ti{s}+dist.overPU(2:end,s))',nrSteps*2,1),...
                    'Color',colorblind('r'));
            end
        end
        % debug ---

        if options.verbose
            disp("...check distances");
        end
        % quick check: distance of H(tk) + underPU(tk) is not ok
        if any(any(dist.affine_tp + dist.underPU > 0))
            % provably falsified!
            savedata.tComp = toc(onlyLoop);

            if options.verbose
                fprintf('...falsification detected!\n');
            end

            % debug ---
            if isplot
                for s=1:nrSpecs
                    subplot(1,nrSpecs,s); hold on;
                    scatter(timeStep*(0:nrSteps),...
                        ds(s)+dist.affine_tp(:,s)+dist.underPU(:,s),...
                        blobsize,'g','filled');
                end
            end
            % debug ---

            % check which spec causes falsification
            for s=1:nrSpecs
                if any(dist.affine_tp(:,s) + dist.underPU(:,s) > 0)
                    falsIdx = s; break
                end
            end
            % time where falsification occurs
            tIdx = find( dist.affine_tp(:,falsIdx) ...
                + dist.underPU(:,falsIdx) > 0,1,'first');
            fals.tFinal = timeStep * (tIdx-1);

            % starting point x0 which yields the falsifying trajectory
            fals.x0 = c_X0 + G_X0*sign(l_prop_T{falsIdx}(tIdx,:)*G_X0)';

            % input contribution: different length depending on whether we
            % have feedthrough or not
            u_tIdx = tIdx - 1;
            if any(any(obj.D))
                u_tIdx = tIdx;
            end
            % contribution from input set
            fals.u = fliplr(generators(origInput.U) * ...
                sign(l_prop_T{falsIdx}(1:u_tIdx,:)*underPU_G)');
            if isuconst
                % contribution from constant input vector
                fals.u = fals.u + origInput.u;
            else
                % contribution from piecewise-constant input vectors
                for j=1:u_tIdx
                    idx = find((j-1)*timeStep >= options.tu,1,'last');
                    fals.u(:,j) = fals.u(:,j) + origInput.u(:,idx);
                end
            end
            if tIdx == 1
                fals.tu = 0;
            else
                fals.tu = (0:timeStep:(tIdx-2)*timeStep)';
            end
            % save data (plots)
            savedata.fals_tFinal = fals.tFinal;
            savedata.dist_affinetp = dist.affine_tp;
            savedata.dist_affinetp_underPU = dist.affine_tp + dist.underPU;
            return;

        else

            % loop over specs...
            timeStep_prop = Inf(nrSteps,nrSpecs);
            for s=1:nrSpecs

                % quick check: all distances H(tauk) + overPU(tk+1) are ok
                if all(dist.affine_ti{s} + dist.overPU(2:end,s) < 0)
                    specUnsat{s} = double.empty(2,0); continue
                end
                
                % check distances one-by-one
                t = [0,timeStep];
                for k=1:nrSteps
        
                    % 1. distance of H(tauk) + overPU(tk+1) is ok
                    if all(dist.affine_ti{s}(k) + dist.overPU(k+1,s) < 0)
                        % add to already verified time intervals
                        specUnsat{s} = aux_removeFromUnsat(specUnsat{s},t);
                    
                    % 2. distance of H(tauk) + overPU(tk+1) is not ok
                    else % case: dist.affine_ti(k) + dist.overPU(k+1) >= 0
                        % suggest refinement for Delta t by using quadratic
                        % dependence of the size of overPU (and Cbloat...)
%                         timeStep_prop(k,s) = aux_affinePUTimeStep(obj,...
%                             dist,options.tFinal,timeStep,timeStepFactor,...
%                             isuconst,options.tu,options.normC,k);
                    end
        
                    % shift time interval
                    t = t + timeStep;
                end

            end

            % method from paper: fixed factor
            timeStep = timeStep * timeStepFactor_fixed;
            % refine time step size
%             timeStep = min(min(timeStep_prop));

            % TODO: shorten time horizon, but ensure that timeStep still
            % yields an integer number of steps

        end

    end


    % specification is satisfied over entire time horizon
    if all(cellfun(@isempty,specUnsat,'UniformOutput',true))
        break
    % for safety reasons, stop in case time step size becomes too small
    % -> return that no decision (-1) could be obtained!
    elseif timeStep < 1e-12
        res = -1; return
    end
    

end

% verification successful
res = true;

% savedata for visualization
savedata.dist_affinetp = dist.affine_tp;
savedata.dist_affineti = dist.affine_ti;
if isU
    savedata.dist_affinetp_underPU = zeros(nrSteps+1,nrSpecs);
    savedata.dist_affineti_overPU = cell(nrSpecs,1);
    for s=1:nrSpecs
        savedata.dist_affinetp_underPU(:,s) = dist.affine_tp(:,s) + dist.underPU(:,s);
        savedata.dist_affineti_overPU{s,1} = dist.affine_ti{s} + dist.overPU(2:end,s);
    end
end

% measure elapsed time
savedata.tComp = toc(onlyLoop);

end


% Auxiliary functions -----------------------------------------------------

function [G,F] = aux_getSetsFromSpec(spec)
% extract safe sets G and unsafe sets F from the specifications
% (to be deleted (once integration in verify is complete))

G = {}; F = {};

for i = 1:length(spec)
    if strcmp(spec(i).type,'safeSet')
        G{end+1}.set = normalizeConstraints(polytope(spec(i).set),'A');
        G{end}.time = spec(i).time;
    elseif strcmp(spec(i).type,'unsafeSet')
        tmp = normalizeConstraints(polytope(spec(i).set),'A');
        if size(tmp.A,1) > 1
            F{end+1}.set = tmp;
            F{end}.time = spec(i).time;
            if size(F{end}.set.A,1) > size(F{end}.set.A,2)
               F{end}.int = interval(F{end}.set);
               F{end}.isBounded = ~(any(isinf(infimum(F{end}.int))) ...
                                   | any(isinf(supremum(F{end}.int))));
            else
               F{end}.isBounded = false; 
            end
        else
            G{end+1}.set = polytope(-tmp.A,-tmp.b);
            G{end}.time = spec(i).time;
        end
    else
       error('This type of specification is not supported!'); 
    end
end

end


% initializations
function [options,origInput] = aux_canonicalForm(obj,options)
% put inhomogeneity to canonical forms:
%    Ax + Bu + c + w  ->  Ax + u, where u \in U + uTransVec
%    Cx + Du + k + v  ->  Cx + v, where v \in V + vTransVec
% the sets options.U and options.V return being centered at the origin, all
% (potentially piecewise-constant) offsets are comprised in the vectors
% options.uTransVec and options.vTransVec

if isa(options.W,'interval')
    options.W = zonotope(options.W);
end
if isa(options.V,'interval')
    options.V = zonotope(options.V);
end

% read out disturbance
centerW = center(options.W);
W = options.W + (-centerW);
% read out sensor noise, combine with feedthrough if given
if any(any(obj.D))
    options.V = obj.D * options.U + options.V;
end
centerV = center(options.V);
options.V = options.V + (-centerV);

% initialize input vector for state and output equation (if sequence given)
if isfield(options,'uTransVec')
    % time-varying input vector
    uVec = options.uTransVec;
else
    % no time-varying input vector, but uTrans given
    uVec = options.uTrans;
end

% save original input for falsifying trajectory
origInput.U = options.U;
origInput.u = uVec;

% put output equation in canonical form
if any(any(obj.D))
    options.vTransVec = obj.D * uVec + obj.k + centerV;
else
    options.vTransVec = obj.k + centerV;
end

% put state equation in canonical form
options.U = obj.B * options.U + W;
options.uTransVec = obj.B * uVec + obj.c + centerW;

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
% initialize one-step propagation matrix
expmat.Deltatk = [];

% precompute inverse of A matrix
expmat.isAinv = rank(full(obj.A)) == obj.dim;
expmat.Ainv = [];
if expmat.isAinv
    expmat.Ainv = inv(obj.A);
end

end

function [dist,fals,isU,G_U,isu,isuconst,options] = aux_initStructsFlags(obj,options)

% struct for distances
dist = [];

% saving of operations for affine systems (u = const. over entire time
% horizon) vs. system with varying u or even uncertainty U
% -> the resulting if-else branching looks quite ugly, but still yields
% large speed-ups for high-dimensional systems
G_U = generators(options.U);
isU = ~isempty(G_U);
isu = any(any(options.uTransVec));
isuconst = size(options.uTransVec,2) == 1;
% sparsity for speed up (acc. to numeric tests only for very sparse
% matrices actually effective)
if nnz(obj.A) / numel(obj.A) < 0.1
    obj.A = sparse(obj.A);
end

% struct for falsifying trajectory
fals.tFinal = []; fals.x0 = []; fals.u = []; fals.tu = [];

end

function [uTrans,vTransNext,tnextSwitch] = aux_uTrans_vTrans(t,...
    timeStep,uTransVec,vTransVec,tu,tnextSwitch,uTrans,vTransNext)
% (note: it is already ensured that the input vector changes over time)
% reads the value of the input vector (uTrans, vTrans) from the matrix
% storing all input vectors (uTransVec, vTransVec) for the current step

% in each step, we compute the output set for
%    R([t_k,t_k+1])   using   vTransVec(k) = vTrans
% and
%    R(t_k+1)   using   vTransVec(k+1) = vTransNext,
% so we have return two values for v

if withinTol(t,tnextSwitch)
    % take correct matrix index depending on current time (only compute
    % this if the next switching time has actually come)
    idx = find(t >= tu,1,'last');
    uTrans = uTransVec(:,idx);

    % update time of next input vector switch
    tnextSwitch = Inf;
    if idx < length(tu)
        tnextSwitch = tu(idx+1);
    end

end
    
if t == 0 || (size(vTransVec,2) > 1 && withinTol(t+timeStep,tnextSwitch))
    % also entered in initial step (to assign vnext), otherwise this
    % if-condition and the above one should not be entered in the same step

    % value for the end of the step
    idx = find(t + timeStep >= tu,1,'last');
    vTransNext = vTransVec(:,idx);

    % no update for next input vector switch here

end

end


% adaptation of time step size
function timeStep = aux_timeStep(timeStep,tFinal,isuconst,tu)
% we compute the largest time step size so that we can use a constant time
% step size which hits all switching times exactly; the input argument
% value for the time step size is an upper bound of the desired value

if isuconst
    % no switches in constant input vector
    timeStep = tFinal / ceil(tFinal / timeStep);
    return
end

% find largest time step size considering switches in piecewise-constant
% input vector

% duration of each piecewise-constant input vector
if withinTol(tFinal-tu(end),0)
    % systems with inhomogeneity in the output equation
    constInt = diff(tu)';
else
    % systems without inhomogeneity in the output equation
    constInt = [diff(tu)',tFinal - tu(end)];
end
% minimum duration and corresponding number of total steps
timeStep = min([min(constInt),timeStep]);
steps = ceil(tFinal / timeStep);

% max number of steps
maxSteps = 1000000000;

% loop over increasing number of steps
while true
    % resulting time step size
    timeStep = tFinal/steps;
    % check if that time step size divides all durations of
    % piecewise-constant input vectors into integers
    temp = constInt ./ timeStep;
    if all(withinTol(temp,round(temp)))
        % time step found
        break
    end

    % increment number of steps
    steps = steps + 1;

    % stopping condition: no time step found until arbitrary value
    if steps > maxSteps
        error("No duration from " + steps + " to " + maxSteps ...
            + " time steps can divide all individual piecewise-constant " ...
            + "input vector durations into integers.");
    end
end

end

function timeStep = aux_affineTimeStep(obj,dist,tFinal,timeStep,...
    timeStepFactor,isuconst,tu,normC)
% adaptation of the time step size for the next iteration based on the
% distance of the affine dynamics solution and the quadratic dependency of
% said distance on the time step size

% for sanity check below...
timeStep_prev = timeStep;

% one time step proposition for each spec
nrSpecs = length(dist.affine_ti);
% alternative using 3 pages
% nrSpecs = size(dist.affine_ti,3);
timeStep_prop = zeros(1,nrSpecs);

for s=1:nrSpecs
    % get time index + start/end of maximum infringement
    [tempDist,tempIdx] = max(dist.affine_ti{s},[],2);
    % alternative using 3 pages
%     [tempDist,tempIdx] = max(dist.affine_ti{s}(:,:,s),[],2);
    [~,tIdx] = max(tempDist);
    
    % compute coefficient for equation: sizeC = (det(eAt)factor)*a*timeStep^2
    % since error set Cbloat shrinks quadratically over timeStep
    sizeFactor = normC * exp(trace(obj.A * timeStep*(tIdx-1)));
    if sizeFactor > 0
        coeff_a = dist.Cbloat(tIdx,s) / (sizeFactor * timeStep^2);
    else
        coeff_a = dist.Cbloat(tIdx,s) / timeStep^2;
    end
    
    % compute time step size which would achieve the maximum
    % admissible value for dist.Cbloat(tIdx)
    tIdx = tIdx + (tempIdx(tIdx) - 1);
    timeStep = sqrt(-dist.affine_tp(tIdx,s) / coeff_a);
    % sanity check
    assert(timeStep < timeStep_prev,'Error in adaptation of time step size');
    
    % adjust time step size to yield an integer number of steps and is
    % compatible with switches in the input vector
    timeStep_prop(s) = aux_timeStep(timeStep,tFinal,isuconst,tu);
end

% ensure that too small values are avoided (which mostly occurs when the
% quadratic approximation does not provide a good estimate, which in turn
% often occurs when the previous time step size was too large)
timeStep_min = aux_timeStep(timeStep_prev*timeStepFactor,tFinal,isuconst,tu);

% take maximum of both time step sizes
timeStep = max([timeStep_prop, timeStep_min]);

end

function timeStep = aux_affinePUTimeStep(obj,dist,tFinal,timeStep,...
    timeStepFactor,isuconst,tu,normC,k)
% adaptation of the time step size for the next iteration based on the
% distance of the full solution and the quadratic dependency of both parts
% of said distance (affine + particular solutions) on the time step size;
% both contributions are affected simultaneously by the same time step size

% for sanity check below...
timeStep_prev = timeStep;

% case: dist.affine_ti(k) + dist.overPU(k+1) >= 0
nrSpecs = length(dist.affine_ti);
timeStep_prop = zeros(1,nrSpecs);

for s=1:nrSpecs

    % compute margin
    surplus = max(dist.affine_tp(k:k+1,s)) + dist.underPU(k+1,s);
    
    % scaling factor by propagation matrix
    sizeFactor = normC * exp(trace(obj.A * timeStep*(k-1)));
    
    % compute coefficient for regression equations:
    % 1. sizeC = (det(eAt)factor)*a*timeStep^2
    % 2. sizePU = (det(eAt)factor)*a*timeStep^2
    % since Cbloat and PU shrink quadratically with the time step size
    if withinTol(sizeFactor,0)
        coeff_a_affine = dist.Cbloat(k,s) / timeStep^2;
        coeff_a_overPU = dist.overPU(k+1,s) / timeStep^2;
    else
        coeff_a_affine = dist.Cbloat(k,s) / (sizeFactor * timeStep^2);
        coeff_a_overPU = dist.overPU(k+1,s) / (sizeFactor * timeStep^2);
    end
    
    % compute time step size which would achieve the maximum
    % admissible value for dist.Cbloat(tIdx,s)
    timeStep_prop(s) = sqrt(-surplus / (coeff_a_affine + coeff_a_overPU));
    
    % adjust time step size to yield an integer number of steps
    timeStep_prop(s) = aux_timeStep(timeStep_prop(s),tFinal,isuconst,tu);

end

% ensure that too small values are avoided (which mostly occurs when the
% quadratic approximation does not provide a good estimate, which in turn
% often occurs when the previous time step size was too large)
timeStep_min = aux_timeStep(timeStep_prev*timeStepFactor,tFinal,isuconst,tu);

% take maximum of both time step sizes
timeStep = max([timeStep_prop, timeStep_min]);

% sanity check
assert(timeStep < timeStep_prev,'Error in adaptation of time step size');

end


% specification-related functions
function FGunsat = aux_removeFromUnsat(FGunsat,t)
% adapt FGunsat so that timeInterval is not part of time intervals covered

FGunsat_col = reshape(FGunsat',numel(FGunsat),1);
if mod(sum(t(1) >= FGunsat_col),2) ~= 0
    % lower bound starts inside unverified time interval
    % t(1) \in [ timeInterval )
    idx = find(t(1) >= FGunsat(:,1) & t(1) <= FGunsat(:,2));
    
    if t(2) <= FGunsat(idx,2)
        if t(1) > FGunsat(idx,1)
            FGunsat = [FGunsat(1:idx-1,:); ...
                [FGunsat(idx,1), t(1)]; FGunsat(idx:end,:)];
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


% interval matrices
function [expmat,intmatF,intmatG] = aux_intmat(obj,isu,expmat,timeStep)
% computation of F and G for a given time step size; currently, these
% variables are computed by a Taylor series until floating-point precision,
% i.e., we increase the truncation order until the additional values are so
% small that the stored number (finite precision!) does not change anymore

% skip computation of Ftilde if u is all-zero

% load data from object/options structure
A = obj.A;
n = obj.dim;

% initialize auxiliary variables and flags for loop
Asum_pos_F = zeros(n);
Asum_neg_F = zeros(n);
stoploop_F = false;
if isu
    Asum_pos_G = zeros(n);
    Asum_neg_G = zeros(n);
    stoploop_G = false;
else
    stoploop_G = true;
    intmatG.set = interval(zeros(n),zeros(n));
    intmatG.rad = zeros(n);
    intmatG.center = zeros(n);
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
            aux_getAposneg(eta,expmat.Apos,expmat.Aneg,Apower_eta);
        
        % if new term does not change result anymore, loop to be finished
        Asum_add_pos_F = factor*expmat.Aneg{eta};
        Asum_add_neg_F = factor*expmat.Apos{eta};
        
        % safety check (if time step size too large, then the sum converges
        % too late so we already have Inf values)
        % previous criterion: any(any(isinf(Asum_add_pos_F))) || any(any(isinf(Asum_add_neg_F)))
        % ... but very costly for large A!
        if eta == 75
            intmatF.set = []; intmatG.set = [];
            expmat.conv = false; return
%             throw(MException('reach_adaptive:notconverging',...
%                 'Time Step Size too big for computation of F.'));
        end
        
        % compute ratio for floating-point precision
        if all(all(Asum_add_pos_F <= eps * Asum_pos_F)) && ...
                all(all(Asum_add_neg_F >= eps * Asum_neg_F))
            stoploop_F = true;
            intmatF.rad = 0.5*(Asum_pos_F - Asum_neg_F);
            intmatF.center = Asum_neg_F + intmatF.rad;
        end

        % compute powers; factor is always negative
        Asum_pos_F = Asum_pos_F + Asum_add_pos_F; 
        Asum_neg_F = Asum_neg_F + Asum_add_neg_F;

    end
    
    if ~stoploop_G
        [expmat.Apos,expmat.Aneg] = ...
            aux_getAposneg(eta-1,expmat.Apos,expmat.Aneg,expmat.Apower{eta-1});
        
        % if new term does not change result anymore, loop to be finished
        % we require one additional division by eta as the terms in expmat
        % are divided by (eta-1)! instead of eta! as required
        Asum_add_pos_G = factor*expmat.Aneg{eta-1} / eta;
        Asum_add_neg_G = factor*expmat.Apos{eta-1} / eta;
        
        % safety check (if time step size too large, then the sum converges
        % too late so we already have Inf values)
        % previous criterion: any(any(isinf(Asum_add_pos_Ftilde))) || any(any(isinf(Asum_add_neg_Ftilde)))
        % ... but very costly for large A!
        if eta == 75
            intmatF.set = []; intmatG.set = [];
            expmat.conv = false; return
%             throw(MException('reach_adaptive:notconverging',...
%                 'Time Step Size too big for computation of Ftilde.'));
        end
        
        % compute ratio for floating-point precision
        if all(all(Asum_add_pos_G <= eps * Asum_pos_G)) && ...
                all(all(Asum_add_neg_G >= eps * Asum_neg_G)) 
            stoploop_G = true;
            intmatG.rad = 0.5*(Asum_pos_G - Asum_neg_G);
            intmatG.center = Asum_neg_G + intmatG.rad;
        end

        % compute powers; factor is always negative
        Asum_pos_G = Asum_pos_G + Asum_add_pos_G; 
        Asum_neg_G = Asum_neg_G + Asum_add_neg_G;
    end
    
    % instantiate interval matrices if converged
    if stoploop_F
        intmatF.set = interval(Asum_neg_F,Asum_pos_F);
    end
    if stoploop_G && isu
        intmatG.set = interval(Asum_neg_G,Asum_pos_G);
    end

    % exit loop if both converged
    if stoploop_F && stoploop_G
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

function [Apos,Aneg] = aux_getAposneg(eta,Apos,Aneg,Apower_eta)
% the separation of A^eta into positive and negative indices can be
% precomputed and saved; Apower_eta has to match eta correctly!

if length(Apos) >= eta && ~isempty(Apos{eta})
    % ... then also length(Aneg) >= eta
    % read from memory
    
else
    
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
            expmat.conv = false; return
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

function [G_underPU,expmat] = aux_underPU(obj,G_U,expmat,timeStep)
% computation of the particular solution due to the input vector u using
% a Taylor series where the truncation order is increased until the
% additional values are so small that the stored number (finite precision!)
% does not change anymore; in case the inverse of A exists, we directly
% compute the analytical solution (where the exponential matrix is also
% only computed until finite precision)
% note: G_U is non-zero (checked outside)


if expmat.isAinv
    G_underPU = expmat.Ainv * (expmat.Deltatk - eye(obj.dim)) * G_U;
    
else    
    % compute by sum until floating-point precision
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
            expmat.conv = false; return
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
    G_underPU = Asum * G_U;
end

end

function [G_overPU,expmat] = aux_overPU(obj,G_U,expmat,timeStep)
% computation of the particular solution due to the uncertain input set U
% using a Taylor series where the truncation order is increased until
% the additional values are so small that the stored number (finite
% precision!) does not change anymore;
% we use only the generator matrix to spare some zonotope function calls

% maximum truncation order
maxeta = 75;

% initialize particular solution
G_U_size = size(G_U,2);
G_overPU = zeros(obj.dim,G_U_size*maxeta);
G_overPU(:,1:G_U_size) = timeStep * G_U;
PU_diag = sum(abs(G_overPU),2);

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
    
    % check if floating-point precision reached
    if all( abs(addG_PU_diag) <= eps * abs(PU_diag) )
        stoploop = true;
    end
    
    % add term to simplified value for convergence
    PU_diag = PU_diag + addG_PU_diag;
    
    % append to generator matrix
    G_overPU(:,eta*G_U_size+1:(eta+1)*G_U_size) = addG_PU;
    
    % break loop, remove remaining pre-allocated zeros
    if stoploop
        G_overPU(:,(eta+1)*G_U_size+1:end) = [];
        break;
    end
    
    % increment eta
    eta = eta + 1;

    % exit and try again with smaller time step size if truncation order
    % becomes too large
    if eta == maxeta
        expmat.conv = false; return
    end

end

end


% ------------------------------ END OF CODE ------------------------------
