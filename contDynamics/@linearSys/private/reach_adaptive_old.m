function [Rout,Rout_tp,res,tVec] = reach_adaptive_old(obj,options)
% reach_adaptive_old - computes the reachable set for linear systems using
%    the propagation from the start as well as adaptive parametrization
%    note: old version (just for bug fixing)
%
% Syntax:
%    [Rout,Rout_tp,res,tVec] = reach_adaptive_old(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    Rout - output set of time intervals for the continuous dynamics 
%    Rout_tp - output set of points in time for the continuous dynamics
%    res - boolean (only if specification given)
%    tVec - vector containing time step sizes (e.g. for plotting)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Mark Wetzlinger
% Written:       23-July-2019
% Last update:   08-Oct-2019
% Last revision: ---


%------------- BEGIN CODE --------------

% setup for calculation of Rout(_tp) --------------------------------------
[C,D,k] = initOutputEquation(obj,options);
% -------------------------------------------------------------------------

% initialize parameters ---------------------------------------------------
options = rA_initTimeStepTaylorTerms(obj,options);
options.zonotopeOrder = Inf;
formerTSTT = [options.timeStep; options.taylorTerms];
tVec = options.timeStep;
% -------------------------------------------------------------------------

% init step ---------------------------------------------------------------
[R_0, options] = rA_initReach(obj, options);
% compute output set
[Rout{1},Rout_tp{1}] = outputSet(C,D,k,R_0,options);
% safety property check
if isfield(options,'specification')
    if ~check(options.specification,Rout{1})
        % violation
        res = false;
        return
    end
end
% -------------------------------------------------------------------------

% access to values for loop -----------------------------------------------
Rhom_0 = options.Rhom;
Rhom_tp_0 = options.Rhom_tp;
eADelta = options.Q;
% -------------------------------------------------------------------------

% main loop: post ---------------------------------------------------------
cycle = 20;
i = 1; % init iteration counter
while options.tFinal - options.t > eps
    % update iteration counter
    i = i + 1;
    if isfield(options,'verbose') && options.verbose
        if mod(i,cycle) == 0
            fprintf("Step " + i + ": " + options.t + "\n");
        end
    end
    
    % find new time step / Taylor terms -----------------------------------
    % updated: options.timeStep|taylorTerms|t
    %          obj.taylor.powers|E|F
    options = rA_postTimeStepTaylorTerms(obj, options);
    tVec(i,1) = options.timeStep;
    % ---------------------------------------------------------------------
    
    
    % calculate new init procedure for current step -----------------------
    if any([options.timeStep; options.taylorTerms] ~= formerTSTT)
        [options, eADelta] = rA_post(obj, options);
        % access values for propagation
        Rhom_0 = options.Rhom;
        Rhom_tp_0 = options.Rhom_tp;
        formerTSTT = [options.timeStep; options.taylorTerms];
    end
    % ---------------------------------------------------------------------
    
    
    % propagate homogeneous part ------------------------------------------
    % H^R([t_i, t_i + Delta t_i]) = e^At_i * H^R([0, Delta t_i])
    Rhom    = options.Q * Rhom_0;
    % H^R([0, Delta t_i) = e^At_i * H^R(Delta t_i)
    Rhom_tp = options.Q * Rhom_tp_0;
    % ---------------------------------------------------------------------
    
    % propagate inhomogeneous part ----------------------------------------
    options = rA_propInhom(obj, options);
    % ---------------------------------------------------------------------
    
    % calculate full reachable set ----------------------------------------
    % R([t_i, t_i + Delta t_i]) = H([t_i, t_i + Delta t_i]) + P([0, t_i])
    Rcont.ti = Rhom + options.Rinhom;
    % R(t_i + Delta t_i) = H(t_i + Delta t_i) + P([0, t_i])
    Rcont.tp = Rhom_tp + options.Rinhom;
    % ---------------------------------------------------------------------
    
    
    % output set and safety property check --------------------------------
    [Rout{i,1},Rout_tp{i,1}] = outputSet(C,D,k,Rcont,options);
    % safety property check
    if isfield(options,'specification')
        if ~check(options.specification,Rout{i})
            % violation
            Rout = Rout(1:i);
            Rout_tp = Rout_tp(1:i);
            res = false;
            return
        end
    end
    % ---------------------------------------------------------------------
    
    
    % propagate matrix exponentials ---------------------------------------
    options.P = options.Q;
    options.Q = options.Q * eADelta;
    % ---------------------------------------------------------------------
    
end


% specification fulfilled at all times
res = true;

end


% Auxiliary Functions -----------------------------------------------------
function options = rA_postTSTT_savedSets(obj, options, idxTimeStep)
% rA_postTSTT_savedSets - calculates
%    obj.taylor.error|powers|F|(inputF)|inputCorr
%    and saves them to
%    options.savedSets{idxTimeStep}.<...>
% note: options.factor calculated in rA_postTSTT_checkTimeStep
%
% Syntax:  
%    options = rA_postTSTT_savedSets(obj, options, idxTimeStep)
%
% Inputs:
%    obj         - linearSys object
%    options     - options struct
%    idxTimeStep - index of current time step in options.savedSets
%
% Outputs:
%    options - options struct
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      31-Aug-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% flag pair of time step and Taylor terms
options.savedSets{idxTimeStep}.eta(options.taylorTerms,1) = true;


% compute error in exponential matrix -------------------------------------
% returns: obj.taylor.powers|error
obj = exponential(obj, options);
options.savedSets{idxTimeStep}.powers = obj.taylor.powers;
options.savedSets{idxTimeStep}.E{options.taylorTerms,1} = obj.taylor.error;
if options.isInput
    options.savedSets{idxTimeStep}.Etu{options.taylorTerms,1} = ...
        obj.taylor.error * options.timeStep * options.u;
end
% -------------------------------------------------------------------------


% compute time interval error (tie) ---------------------------------------
% returns: obj.taylor.F
obj = tie(obj,options);
options.savedSets{idxTimeStep}.F{options.taylorTerms,1} = obj.taylor.F;
options.savedSets{idxTimeStep}.FR0{options.taylorTerms,1} = ...
    obj.taylor.F * options.R0;
% -------------------------------------------------------------------------


% obtain factors for inputTie ---------------------------------------------
%  ...only if 0 not in input set
if ~options.originContained
    % compute reachable set due to input:
    % returns: obj.taylor.inputF|inputCorr
    % ...calc. obj.taylor.V|RV|Rtrans|eAtInt|timeStep|eAt in syn_post!
    obj = inputTie(obj,options);
    if isempty(obj.c); vTrans = obj.B*options.uTrans;
    else;              vTrans = obj.B*options.uTrans + obj.c; end
    obj.taylor.inputCorr = obj.taylor.inputF * zonotope(vTrans);
    options.savedSets{idxTimeStep}.inputCorr{options.taylorTerms,1} = ...
        obj.taylor.inputCorr;
else
    options.savedSets{idxTimeStep}.inputCorr{options.taylorTerms,1} = ...
        zeros(obj.dim,1);
end
% -------------------------------------------------------------------------



end

function [idxTimeStep, options] = rA_postTSTT_checkTimeStep(obj, options)
% rA_postTSTT_checkTimeStep - check if current time step has been used
%    before and returns index for options.savedSets
% furthermore, if the time step is new and 0 \notin input set,
%    options.factor until options.taylorTerms + 1 is initialized
%    (this is done to save operations)
%
% Syntax:  
%    [idxTimeStep, options] = rA_postTSTT_checkTimeStep(options)
%
% Inputs:
%    options - options struct
%
% Outputs:
%    idxTimeStep - index of current time step in options.savedSets
%    options - options struct
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      31-Aug-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

idxTimeStep = find(abs(options.usedTimeSteps - options.timeStep) < 1e-9,1); 
if isempty(idxTimeStep)
    % new time step (find returns [] if not found)
    idxTimeStep = length(options.savedSets) + 1;
    options.usedTimeSteps(idxTimeStep,1) = options.timeStep;
    options.savedSets{idxTimeStep,1}.timeStep = options.timeStep;
    options.savedSets{idxTimeStep,1}.eta = false(options.maxTaylorTerms,1);
    
    % initialize options.factor
    options.factor = 0;
    % used in inputTie and inputSolution
    for i=1:(options.maxTaylorTerms+1)
        % compute initial state factor
        options.factor(i) = options.timeStep^(i)/factorial(i);    
    end
    options.savedSets{idxTimeStep}.factor = options.factor;
end

end

function options = rA_postTimeStepTaylorTerms(obj, options)
% rA_postTimeStepTaylorTerms - find time step and
%  number of Taylor terms for new iteration
%
% Syntax:  
%    options = rA_postTimeStepTaylorTerms(obj, options)
%
% Inputs:
%    obj     - linearSys object
%    options - options struct
%
% Outputs:
%    options - options struct
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      25-Aug-2019
% Last update:  08-Oct-2019
% Last revision:---

%------------- BEGIN CODE --------------


remainingTime = options.tFinal - options.t;
if remainingTime <= options.timeStep * options.fraction
    % use remaining time
    options.timeStep = remainingTime;
else
    % try bigger time step
    options.timeStep = options.timeStep * options.fraction;
end
options.taylorTerms = options.maxTaylorTerms - 1; % reset
maxEnough = true;

[idxTimeStep, options] = rA_postTSTT_checkTimeStep(obj,options);
% initialize bounds -------------------------------------------------------
if options.isInput
    remainingSteps = ceil((options.tFinal - options.t)/options.timeStep);
    % allowed error for this step: allowed - curr divided by forseen time
    options.ePadm = (options.ePmax - options.eP) / remainingSteps;
    % ... if no inhom. solution, then ePadm same as in init
end
% -------------------------------------------------------------------------

% loop until bloatings due to Q * F * R0 and Q * E * t * u ...
%    below their respective allowed thresholds ----------------------------
while true
    
    options.taylorTerms = options.taylorTerms + 1;
    if options.taylorTerms > options.maxTaylorTerms
        if maxEnough
            options.taylorTerms = 0;
            maxEnough = false;
            continue
        end
        % if max eta reached, shrink time step and try again
        options.timeStep = options.timeStep / options.fraction;
        % recompute bound for QEtu
        if options.isInput
            remainingSteps = ceil((options.tFinal - options.t)/options.timeStep);
            % allowed error for this step: allowed - curr divided by forseen time
            options.ePadm = (options.ePmax - options.eP) / remainingSteps;
            % ... if no particular solution, then ePadm same as in init
        end
        [idxTimeStep, options] = rA_postTSTT_checkTimeStep(obj,options);
        options.taylorTerms = options.maxTaylorTerms - 1;
        maxEnough = true; % hopeful set for smaller time step
        continue % start with new time step and eta = max new while iteration
    end
    
    if ~options.savedSets{idxTimeStep}.eta(options.taylorTerms)
        % calculates obj.taylor.error|powers|F|inputCorr and options.factor
        % .. and saves them to options.savedSets{idxTimeStep}
        options = rA_postTSTT_savedSets(obj, options, idxTimeStep);
    end
    
    % calculate hom. error
    if ~options.originContained
        eH = norm(sum(abs(generators(options.Q * ...
                options.savedSets{idxTimeStep}.FR0{options.taylorTerms})),2),2) + ...
            norm(sum(abs(generators(options.Q * ...
                options.savedSets{idxTimeStep}.inputCorr{options.taylorTerms})),2),2);
    else
        eH = norm(sum(abs(generators(options.Q * ...
                options.savedSets{idxTimeStep}.FR0{options.taylorTerms})),2),2);
    end
    
    if eH > options.eHmax
        % max Taylor terms does not satisfy hom. bloating
        maxEnough = false;
        continue
    end
    
    if options.isInput
        addErrorP = norm(sum(abs(generators(options.Q * ...
            options.savedSets{idxTimeStep}.Etu{options.taylorTerms})),2),2);
    else
        addErrorP = norm(sum(abs(generators(options.Q * ...
            (options.savedSets{idxTimeStep}.E{options.taylorTerms} * options.R0))),2),2);
    end
    
    if addErrorP > options.ePadm
        % max Taylor terms does not satisfy hom. bloating
        maxEnough = false;
        continue
    end
    
    % if both below and not maxEnough
    if ~maxEnough && eH < options.eHmax && addErrorP < options.ePadm
         break;
    end

end % loop exit if both timeStep and Taylor terms are ok for F and Inhom bloating
% -------------------------------------------------------------------------

% insert correct values in obj and options struct -------------------------
if options.isInput
    options.Etu = options.savedSets{idxTimeStep}.Etu{options.taylorTerms};
end
options.FR0 = options.savedSets{idxTimeStep}.FR0{options.taylorTerms};

% for inputSolution
if ~options.originContained
    options.factor = options.savedSets{idxTimeStep}.factor;
    obj.taylor.inputCorr = options.savedSets{idxTimeStep}.inputCorr{options.taylorTerms};
end

obj.taylor.powers = options.savedSets{idxTimeStep}.powers;
obj.taylor.error = options.savedSets{idxTimeStep}.E{options.taylorTerms};
obj.taylor.F = options.savedSets{idxTimeStep}.F{options.taylorTerms};
% obj.taylor is correct with the exception of:
%    obj.taylor.V|RV|Rtrans|eAtInt|timeStep|eAt
% ...currently, obj.taylor.inputF|inputCorr are also calculated in post
% -------------------------------------------------------------------------

% accumulate bloating -----------------------------------------------------
if options.isInput
    addErrorP = sum(abs(generators(options.Q * options.Etu)),2);
    options.ePbox = options.ePbox + addErrorP;
    options.eP = norm(options.ePbox,2);
end
% -------------------------------------------------------------------------

% propagate current time --------------------------------------------------
options.t = options.t + options.timeStep;
% -------------------------------------------------------------------------

end

function [options, eAt] = rA_post(obj, options)
% rA_post - post step for adaptively parametrized reachability analysis
%
% Syntax: [options, eAt] = rA_post(obj, options)
%
% Inputs:
%    obj      - linearSys object
%    options  - options struct, containing
%                    .R0: initial set
%                    .timeStep: time step size for current iteration
%                    .taylorTerms: taylor Terms for current iteration
%
% Outputs:
%    options - options struct, containing
%                    .Rhom: homogeneous time interval solution
%                    .Rhom_tp: homogeneous time point solution
%                    .Rtrans: inhomogeneous solution due to uTrans
%                    .Raux: inhomogeneous solution due to U
%    eAt     - propagation matrix to end of current time step
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      25-Aug-2019
% Last update:  31-Aug-2019
% Last revision:---

%------------- BEGIN CODE --------------

% preliminary computations ------------------------------------------------
% compute exponential matrix:
% returns: obj.taylor.powers|E
% obj = exponential(obj,options);

% compute time interval error (tie):
% returns: obj.taylor.F
% obj = tie(obj,options);

% ^ both above already calculated in rA_init_TimeStepTaylorTerms

% compute reachable set due to input:
% assert correct obj.taylor.error|powers
% returns: obj.taylor.V|RV|Rtrans|eAtInt
% ...also obj.taylor.inputF|inputCorr (already calculated in postTSTT...)
obj = inputSolution(obj,options);

% save time step using for above calculations in same struct
obj.taylor.timeStep = options.timeStep;

% compute reachable set of first time interval
eAt = expm(obj.A*options.timeStep);

% save current propagation matrix in same struct
obj.taylor.eAt = eAt;
% -------------------------------------------------------------------------


% read necessary sets from preliminary computations -----------------------
% assert obj.taylor.F|RV|Rtrans|inputCorr
F = obj.taylor.F;
if any(any(options.U.Z))
    gensRV = generators(obj.taylor.RV);
    RV = zonotope([zeros(obj.dim,1), gensRV(:,any(gensRV,1))]);
else
    RV = zonotope(zeros(obj.dim,1));
end
if any(options.uTrans)
    if iscell(obj.taylor.Rtrans)
        Rtrans = obj.taylor.Rtrans{1};
    else
        gensRtrans = generators(obj.taylor.Rtrans);
        Rtrans = zonotope([center(obj.taylor.Rtrans), ...
            gensRtrans(:,any(gensRtrans,1))]);
    end
else
    Rtrans = zonotope(zeros(obj.dim,1));
end
inputCorr = obj.taylor.inputCorr;
Rinit = options.R0;
% -------------------------------------------------------------------------


% first time step homogeneous solution ------------------------------------
Rhom_tp = eAt*Rinit + Rtrans;
if isa(Rinit,'quadZonotope')
    Rhom = enclose(Rinit,Rhom_tp) + F*zonotope(Rinit) + inputCorr;
elseif isa(Rinit,'zonoBundle') 
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit.Z{1} + inputCorr;
else
    enc = enclose(Rinit,Rhom_tp);
    bloat = F*Rinit;
    Rhom = enc + bloat + inputCorr;
end
% note: no reduction
% -------------------------------------------------------------------------

% save homogeneous and particulate solution -------------------------------
options.Rhom    = Rhom;
options.Rhom_tp = Rhom_tp;
options.Raux    = RV;
options.Rtrans  = Rtrans;
% -------------------------------------------------------------------------

end

function options = rA_propInhom(obj, options)
% rA_propInhom - propagate inhomogeneous solution,
%    including reduction of zonotope order until limit given by
%    allowed bloating reached
% formula given by:
% P^R([0, t_i + Delta t_i]) = P^R([0, t_i])
%    + e^At_i * P^R_U(Delta t_i) + e^At_i-1 * P^R_u(Delta t_i-1)
%
% Syntax:  
%    options = rA_propInhom(obj,options)
%
% Inputs:
%    obj     - linearSys object
%    options - options struct
%
% Outputs:
%    options - options struct
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      29-Aug-2019
% Last update:  15-October-2019
% Last revision:---

%------------- BEGIN CODE --------------

% calculate linearly growing boundary for allowed bloating error
options.eSadm = options.eSmax * options.t / options.tFinal;

% additional inhomogenuity due to current time step
Rinhomadd = options.Q * options.Raux + options.P * options.lastRtrans;

if options.isInput
    % generators (+length) of added part
    options.RinhomG = [options.RinhomG, generators(Rinhomadd)];
    options.RinhomGlength = [options.RinhomGlength, vecnorm(generators(Rinhomadd),2)];

    % sort generators by length
    numGensInhom = length(options.RinhomGlength);
    [~, idxRinhomG] = mink(options.RinhomGlength,numGensInhom);
    
    binary = true; % binary search vs. start-to-end search
    
    if binary
    
    lower = 0;
    upper = numGensInhom + 1; 
    gen = floor((upper - lower) * 0.5);
    while true
        addErrorS = vecnorm(sum(abs(options.RinhomG(:,idxRinhomG(1:gen))),2),2);
        if options.eS + addErrorS > options.eSadm
            % solution < gen
            if gen == 1
                break
            end
            upper = gen;
        else
            % solution > gen
            if gen == lower
                % gen is max. number that can be over-approximated
                break
            end
            lower = gen;
        end
        % update curr
        gen = lower + floor((upper - lower) * 0.5); 
    end
    
    redIdx = idxRinhomG(1:gen);
    nonredIdx = idxRinhomG(gen+1:numGensInhom);
    options.RinhomConvInt = options.RinhomConvInt + ...
        diag(sum(abs(options.RinhomG(:,redIdx)),2));
    options.eS = norm(sum(abs(options.RinhomConvInt),2),2);
    % carry only non-reduced generators
    options.RinhomG = options.RinhomG(:,nonredIdx);
    options.RinhomGlength = options.RinhomGlength(nonredIdx);
    % add centers, take [non-reduced generators, reduced generators]
    options.Rinhom = ...
        zonotope([center(options.Rinhom)+center(Rinhomadd),...
        options.RinhomG,options.RinhomConvInt]);
    else
    
    for gen=1:numGensInhom
        addErrorS = norm(sum(abs(options.RinhomG(:,idxRinhomG(1:gen))),2),2);
        if options.eS + addErrorS > options.eSadm || gen == numGensInhom
            nonredIdx = idxRinhomG(gen:end);
            options.RinhomConvInt = options.RinhomConvInt + ...
                diag(sum(abs(options.RinhomG(:,idxRinhomG(1:gen-1))),2));
            options.eS = norm(sum(abs(options.RinhomConvInt),2),2);
            % carry only non-reduced generators
            options.RinhomG = options.RinhomG(:,nonredIdx);
            options.RinhomGlength = options.RinhomGlength(nonredIdx);
            % add centers, take [non-reduced generators, reduced generators]
            options.Rinhom = ...
                zonotope([center(options.Rinhom)+center(Rinhomadd),...
                options.RinhomG,options.RinhomConvInt]);
            break
        end
    end

    end
    
    % monitor how many generators are reduced
    options.RinhomRed(end+1,1) = gen;
end

% update Rtrans
options.lastRtrans = options.Rtrans;

end

function options = rA_initTimeStepTaylorTerms(obj, options)
% rA_initTimeStepTaylorTerms - find time step size and taylor terms
%
% Syntax:  
%    options = rA_initTimeStepTaylorTerms(obj, options)
%
% Inputs:
%    obj     - linearSys
%    options - options for the computation of reachable sets
%
% Outputs:
%    options - options for the computation of reachable sets
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      25-Aug-2019
% Last update:  31-Aug-2019
% Last revision:---

%------------- BEGIN CODE --------------


% define shrinking factor -------------------------------------------------
shrinkT = 0.8; % (0, 1)
% -------------------------------------------------------------------------

% initialize values for loop ----------------------------------------------
options.timeStep = (options.tFinal - options.tStart) * 0.1;
% all time steps should not have more Taylor terms than below
options.maxTaylorTerms = 10;
options.taylorTerms = options.maxTaylorTerms - 1;
% initialize options.factor until options.maxTaylorTerms
options.factor = 0;
for i=1:(options.maxTaylorTerms+1)
    %compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end

options.u = obj.B * (options.U + options.uTrans);
options.isInput = any(any(options.u.Z));
if ~options.isInput
    % outside of loop since independent of options.timeStep
    options.eHmax = options.error;
    options.ePmax = 0.01 * options.eHmax;
    options.ePadm = options.ePmax;
    options.eSmax = options.ePmax;
    % init zonotope order error for later
    options.eSadm = options.eSmax;
    options.eS = 0;
else
    options.eHmax = options.error * 0.30;
    options.ePmax = options.error * 0.30;
    remainingSteps = ceil(options.tFinal / options.timeStep);
    options.ePadm = options.ePmax / remainingSteps;
    % init zonotope order error for later
    options.eSmax = options.error * 0.40;
    options.eSadm = options.eSmax / remainingSteps;
    options.eS = 0;
end
% -------------------------------------------------------------------------


% loop until bloatings due to Q * F * R0 and Q * E * t * u ...
%    below their respective allowed thresholds ----------------------------
maxEnough = true;
while true
    
    options.taylorTerms = options.taylorTerms + 1;
    if options.taylorTerms > options.maxTaylorTerms
        if maxEnough
            options.taylorTerms = 0;
            maxEnough = false;
            continue
        end
        % if max eta reached, shrink time step and try again
        options.timeStep = options.timeStep * shrinkT;
        % recompute bound for QEtu
        if options.isInput
            remainingSteps = ceil(options.tFinal / options.timeStep);
            options.ePadm = options.ePmax / remainingSteps;
            options.eSadm = options.eSmax / remainingSteps;
        end
        % initialize options.factor until options.maxTaylorTerms
        options.factor = 0;
        for i=1:(options.maxTaylorTerms+1)
            %compute initial state factor
            options.factor(i) = options.timeStep^(i)/factorial(i);    
        end
        options.taylorTerms = options.maxTaylorTerms - 1;
        maxEnough = true; % init for smaller time step
        continue % start with new time step and eta = max new while iteration
    end
    
    % returns obj.taylor.error|powers
    obj = exponential(obj,options);
    if options.isInput
        options.Etu = obj.taylor.error * options.timeStep * options.u;
        options.ePbox = sum(abs(generators(options.Etu)),2);
        options.eP = norm(options.ePbox,2);
    else
        options.ePbox = 0;
        options.eP = norm(sum(abs(generators(obj.taylor.error * options.R0)),2),2);
    end

    % if inhom error already to big, start loop anew
    if options.eP > options.ePadm
        % max Taylor terms does not satisfy hom. bloating
        maxEnough = false;
        continue
    end
    
    % compute obj.taylor.F ------------------------------------------------
    obj = tie(obj,options);
    options.FR0 = obj.taylor.F * options.R0;
    % ---------------------------------------------------------------------
    
    % obtain factors for inputSolution > inputTie (only if 0 not in input set)
    if ~options.originContained
        % compute reachable set due to input:
        % returns: obj.taylor.inputF|inputCorr
        % ...calc. obj.taylor.V|RV|Rtrans|eAtInt|timeStep|eAt in syn_initReach!
        obj = inputTie(obj,options);
        if isempty(obj.c); vTrans = obj.B*options.uTrans;
        else;              vTrans = obj.B*options.uTrans + obj.c; end
        obj.taylor.inputCorr = obj.taylor.inputF * zonotope(vTrans);
        % compute induced bloating due to F * X0 + Ftilde * zon(uTrans)
        eH = norm(sum(abs(generators(options.FR0)),2),2) + ...
            norm(sum(abs(generators(obj.taylor.inputCorr)),2),2);
    else
        % compute induced bloating due to F * X0 (inputCorr = 0)
        eH = norm(sum(abs(generators(options.FR0)),2),2);
    end
    % ---------------------------------------------------------------------

    % if hom error already to big, start loop anew ------------------------
    if eH > options.eHmax
        maxEnough = false;
        continue
    end
    % ---------------------------------------------------------------------
    
    % if both below and not maxEnough
    if ~maxEnough && eH < options.eHmax && options.eP < options.ePadm
         break;
    end
    
end
% -------------------------------------------------------------------------

% saved sets for look-up instead of recalculation (time save) -------------
options.usedTimeSteps(1) = options.timeStep;
options.savedSets{1}.timeStep = options.timeStep;
options.savedSets{1}.eta = false(options.maxTaylorTerms,1);
options.savedSets{1}.eta(options.taylorTerms) = true;
options.savedSets{1}.powers = obj.taylor.powers;
options.savedSets{1}.E{options.taylorTerms,1} = obj.taylor.error;
if options.isInput
    options.savedSets{1}.Etu{options.taylorTerms,1} = options.Etu;
end
options.savedSets{1}.F{options.taylorTerms,1} = obj.taylor.F;
options.savedSets{1}.FR0{options.taylorTerms,1} = options.FR0;
options.savedSets{1}.factor = options.factor;
if ~options.originContained
    options.savedSets{1}.inputCorr{options.taylorTerms,1} = ...
        obj.taylor.inputCorr;
else
    options.savedSets{1}.inputCorr{options.taylorTerms,1} = ...
        zeros(obj.dim,1);
end
% -------------------------------------------------------------------------

% fill up the rest of the init time step for lower eta --------------------
% save taylor model for syn_initReach
taylorModelFirstStep = obj.taylor;

while options.taylorTerms > 1
    options.taylorTerms = options.taylorTerms - 1;
    
    % caution: below not equivalent to rA_postTSTT_savedSets!
    % ... here: not saving powers|factor
    
    options.savedSets{1}.eta(options.taylorTerms,1) = true;
    
    % returns: obj.taylor.powers|error
    obj = exponential(obj, options);
    % don't override options.savedSets{1}.powers!
    options.savedSets{1}.E{options.taylorTerms,1} = obj.taylor.error;
    if options.isInput
        options.savedSets{1}.Etu{options.taylorTerms,1} = ...
            obj.taylor.error * options.timeStep * options.u;
    end

    % compute time interval error (tie) -----------------------------------
    % returns: obj.taylor.F
    obj = tie(obj,options);
    options.savedSets{1}.F{options.taylorTerms,1} = obj.taylor.F;
    options.savedSets{1}.FR0{options.taylorTerms,1} = ...
        obj.taylor.F * options.R0;
    % ---------------------------------------------------------------------

    % obtain factors for inputTie -----------------------------------------
    %  ...only if 0 not in input set
    if ~options.originContained        
        % compute reachable set due to input:
        % returns: obj.taylor.inputF|inputCorr
        % ...calc. obj.taylor.V|RV|Rtrans|eAtInt in syn_post!
        obj = inputTie(obj,options);
        % compute vTrans 
        vTrans = obj.B*options.uTrans;
        % consider constant input
        if ~isempty(obj.c)
            vTrans = vTrans + obj.c;
        end
        obj.taylor.inputCorr = obj.taylor.inputF * zonotope(vTrans);
        options.savedSets{1}.inputCorr{options.taylorTerms,1} = ...
            obj.taylor.inputCorr;
    else
        options.savedSets{1}.inputCorr{options.taylorTerms,1} = ...
            zeros(obj.dim,1);
    end
end

% reset taylorModel and taylorTerms to correct value
% ... otherwise obj.taylor for options.taylorTerms = 1 used!
obj.taylor = taylorModelFirstStep;
options.taylorTerms = find(options.savedSets{1}.eta,1,'last');
% -------------------------------------------------------------------------


% init current time
options.t = options.timeStep;

% fraction to de-/increase time step size in post:
options.fraction = 1.1; % arbitrarily set

end

function [Rfirst, options] = rA_initReach(obj, options)
% rA_initReach - init step: R([0, Delta t_0])
%
% Syntax:  
%    [Rfirst, options] = rA_initReach(obj, options)
%
% Inputs:
%    obj      - linearSys object
%    options  - options struct, containing
%                    .R0: initial set
%                    .timeStep: time step size for current iteration
%                    .taylorTerms: taylor Terms for current iteration
%
% Outputs:
%    Rfirst  - reachable set for time interval (.ti) and time point (.tp)
%    options - options struct, containing
%                    .Rhom: homogeneous time interval solution
%                    .Rhom_tp: homogeneous time point solution
%                    .isInhom: if inhomogenuity present
%                    .Rtrans: inhomogeneous solution due to uTrans
%                    .Raux: inhomogeneous solution due to U
%                    .Rinhom: full inhomogeneous solution
%                    .P: identity matrix (used for propagation of Rtrans)
%                    .Q: propagation matrix using .timeStep (for Raux)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      25-Aug-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% preliminary computations ------------------------------------------------
% compute exponential matrix:
% returns: obj.taylor.powers|error
% obj = exponential(obj,options);

% compute time interval error (tie):
% returns: obj.taylor.F
% obj = tie(obj,options);

% ^ both above already calculated in init_TimeStepTaylorTerms

if options.isInput
    options.factor = options.savedSets{1}.factor;
end
% compute reachable set due to input:
% returns: obj.taylor.V|RV|Rtrans|inputF|inputCorr|eAtInt
obj = inputSolution(obj,options);

% save time step using for above calculations in same struct
obj.taylor.timeStep = options.timeStep;

% compute reachable set of first time interval
eAt = expm(obj.A*options.timeStep);

% initialize propagation matrices for analysis
options.P = speye(obj.dim);
options.Q = eAt;

% save current propagation matrix in same struct
obj.taylor.eAt = eAt;
% -------------------------------------------------------------------------


% read necessary sets from preliminary computations -----------------------
F = obj.taylor.F;
Rinit = options.R0;
if any(any(options.U.Z))
    gensRV = generators(obj.taylor.RV);
    RV = zonotope([zeros(obj.dim,1), gensRV(:,any(gensRV,1))]);
else
    RV = zonotope(zeros(obj.dim,1));
end
if any(options.uTrans)
    if iscell(obj.taylor.Rtrans)
        Rtrans = obj.taylor.Rtrans{1};
    else
        gensRtrans = generators(obj.taylor.Rtrans);
        Rtrans = zonotope([center(obj.taylor.Rtrans), gensRtrans(:,any(gensRtrans,1))]);
    end
else
    Rtrans = zonotope(zeros(obj.dim,1));
end
inputCorr = obj.taylor.inputCorr;
% -------------------------------------------------------------------------


% first time step homogeneous solution ------------------------------------
Rhom_tp = eAt*Rinit + Rtrans;
if isa(Rinit,'quadZonotope')
    Rhom = enclose(Rinit,Rhom_tp) + F*zonotope(Rinit) + inputCorr;
elseif isa(Rinit,'zonoBundle') 
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit.Z{1} + inputCorr;
else
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit + inputCorr;
end
% note: no reduction
% -------------------------------------------------------------------------


% save homogeneous and inhomogeneous solution -----------------------------
options.Rhom    = Rhom;
options.Rhom_tp = Rhom_tp;

options.Raux    = RV;
options.Rtrans  = Rtrans;
options.lastRtrans = options.Rtrans;

options.Rinhom = RV;
if options.isInput
    % generators of Rinhom
    options.RinhomG = generators(options.Rinhom);
    options.RinhomGlength = vecnorm(options.RinhomG,2);
    options.RinhomRed = 0;
    % propInhom_new ---
%     options.RinhomGorder = vecnorm(options.RinhomG,1) - vecnorm(options.RinhomG,Inf);
    % ---
    % container for generators converted to interval (always diag matrix)
    options.RinhomConvInt = zeros(obj.dim);
end
% -------------------------------------------------------------------------


% total solution ----------------------------------------------------------
if isa(Rinit,'mptPolytope')
    % convert zonotopes to polytopes
    Radd      = mptPolytope(RV);
    Rfirst.ti = Rhom + Radd;
    Rfirst.tp = Rhom_tp + Radd;
else
    % original computation
    Rfirst.ti = Rhom + RV;
    Rfirst.tp = Rhom_tp + RV;
end
% -------------------------------------------------------------------------

end




%------------- END OF CODE --------------