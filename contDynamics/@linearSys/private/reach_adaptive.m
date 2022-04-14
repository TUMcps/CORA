function [Rout,Rout_tp,res,tVec] = reach_adaptive(obj,options)
% reach_adaptive - computes the reachable set for linear systems using the
%  propagation from the start as well as adaptive parametrization
%
% Syntax:
%    [Rout,Rout_tp,res,tVec] = reach_adaptive(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    Rout - output set of time intervals for the continuous dynamics 
%    Rout_tp - output set of points in time for the continuous dynamics
%    res - boolean (only if specification given)
%    tVec - vector containing time step sizes (for reachSet object)
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
% Last revision: 11-February-2021


%------------- BEGIN CODE --------------

% init global iteration counter
options.i = 0;
% set time
options.t = 0;

% setup for calculation of Rout(_tp)
[C,D,k] = initOutputEquation(obj,options);

% special options
options.timeStepScaling = 1.1; % searching for timeStep
options.zetaTlin = 0.0005; % for initTimeStep (per uTrans)

% initialize input vector if sequence given
if isfield(options,'uTransVec')
    [options.uTransVec,options.uTransVecSwitch] = inputTrajectory(obj,options);
    options.uTrans = options.uTransVec(:,1);
    options.isInput = true;
else
    options.uTransVecSwitch = 0;
    % set input options
    fullU = obj.B * (options.U + options.uTrans);
    options.isInput = any(any(center(fullU))) || any(any(generators(fullU)));
end

% init set for loop
Rcont.tp = [];

% update errors
options = errorAllocation(obj,options);

% figure; hold on;

% loop over each entry in uTransVec once with the algorithm
for u=1:length(options.uTransVecSwitch)
    
    % increase iteration counter
    options.i = options.i + 1;
    
    % log information
    verboseLog(options.i,options.t,options);
    
%     plot(options.R0,[1,2],'r');
    
    % update initial set (last time-point solution) and time step
    if ~isempty(Rcont.tp)
        options = updateStartSet(obj,options,Rcont.tp);
    end
    
%     plot(options.R0,[1,2],'r');
    
    % update uTrans
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,u);
    end

    % update tFinal
    if u < length(options.uTransVecSwitch)
        tFinal = options.uTransVecSwitch(u+1);
    else
        tFinal = options.tFinal;
    end
    
    % find init time step size of current uTrans
    options = initTimeStep(obj,options);
    options.zonotopeOrder = Inf;
    tVec(options.i,1) = options.timeStep;

    % init step
    options = rA_initReach_new(obj, options);
    % access to values for loop
    H_0 = options.Rhom;
    H_tp_0 = options.Rhom_tp;
    % compute first reachable set
    R_0.ti = H_0 + options.Rinhom;
    R_0.tp = H_tp_0 + options.Rinhom;
    
%     plot(R_0.tp,[1,2],'b');

    % update time
    options.t = options.t + options.timeStep;
    
    % compute output set
    [Rout{options.i,1},Rout_tp{options.i,1}] = outputSet(C,D,k,R_0,options);

    % safety property check
    if isfield(options,'specification')
        if ~check(options.specification,Rout{1})
            % violation
            res = false;
            return
        end
    end

    % main loop over constant uTrans
    while tFinal - options.t > 1e-9
        % update iteration counter
        options.i = options.i + 1;

        % log information
        verboseLog(options.i,options.t,options);

        % obtain new time step
        options = postTimeStep(obj, options);
        tVec(options.i,1) = options.timeStep;

        % recompute sets if time step size changed
        if options.timeStep ~= tVec(options.i-1)
            % compute new init procedure for current step
            options = rA_initReach_new(obj, options);
            % access values for propagation
            H_0 = options.Rhom;
            H_tp_0 = options.Rhom_tp;
        end

        % propagate homogeneous part
        % H^R([t_i, t_i + Delta t_i]) = e^At_i * H^R([0, Delta t_i])
        H    = options.Q * H_0;
        % H^R([0, Delta t_i]) = e^At_i * H^R(Delta t_i)
        H_tp = options.Q * H_tp_0;

        % propagate inhomogeneous part
        options = propInhom(obj, options);
        % update last Rtrans
        options.lastRtrans = options.Rtrans;

        % calculate full reachable set
        % R([t_i, t_i + Delta t_i]) = H([t_i, t_i + Delta t_i]) + P([0, t_i])
        Rcont.ti = H + options.Rinhom;
        % R(t_i + Delta t_i) = H(t_i + Delta t_i) + P([0, t_i])
        Rcont.tp = H_tp + options.Rinhom;
        
%         plot(Rcont.tp,[1,2],'b');

        % output set
        [Rout{options.i,1},Rout_tp{options.i,1}] = outputSet(C,D,k,Rcont,options);
        
        % safety property check
        if isfield(options,'specification')
            if ~check(options.specification,Rout{options.i})
                % violation
                Rout = Rout(1:options.i);
                Rout_tp = Rout_tp(1:options.i);
                res = false;
                return
            end
        end

        % update time
        options.t = options.t + options.timeStep;

        % propagate matrix exponentials
        options.P = options.Q;
        options.Q = options.Q * obj.taylor.eAt;

    end

end

% log information
verboseLog(options.i+1,options.t,options);

% specification fulfilled at all times
res = true;

end



% Auxiliary Functions -----------------------------------------------------

function [uTransVec,uTransVecSwitch] = inputTrajectory(obj,options)
% find all switching times in uTransVec, rewrite uTransVec if repetitions

% determine switching times
uTransVec = options.uTransVec;
uTransVecSwitch = linspace(options.tStart,options.tFinal,size(uTransVec,2)+1);
uTransVecSwitch = uTransVecSwitch(1:end-1);

% find repeated entries in uTransVec
remIdx = true(1,size(uTransVec,2));
startIdx = 1;
while startIdx < size(uTransVec,2)
    shift = 1;
    while true
        if ~all(abs(uTransVec(:,startIdx) - uTransVec(:,startIdx+shift)) < 1e-9)
            % uTransVec entry at <startIdx+shift> differs from entry <startIdx>
            break;
        end
        remIdx(startIdx+shift) = false;
        shift = shift + 1;
    end
    startIdx = startIdx + shift;
end

% remove repeated entries
uTransVec = obj.B * uTransVec(:,remIdx);
uTransVecSwitch = uTransVecSwitch(:,remIdx);

end

function options = updateStartSet(obj,options,lastRtp)
% updated: options.R0|S

% calculate linearly growing boundary for allowed bloating error
fracRemTime = options.timeStep / (options.tFinal - options.t);
eSadm = (options.eSmax - options.eS) * fracRemTime;

% sort generators by length
G = generators(lastRtp);
numGens = size(G,2);
lengthG = vecnorm(G,2);
[~, idxLength] = mink(lengthG,numGens);

% accumulate errors of uTrans from before
options.ePprev = options.ePprev + options.eP;
options.eSprev = options.eSprev + options.eS;
% reset errors
options.ePbox = 0;
options.eP = 0;
options.eS = 0;

% add error
addErrorSprev = 0;
for gen=1:numGens
    addErrorS = norm(sum(abs(G(:,idxLength(1:gen))),2),2);
    if options.eSprev + addErrorS > eSadm || gen == numGens
        % add error to eSprev, since not same computation as in propInhom
        options.eSprev = options.eSprev + addErrorSprev;
        % add centers, take [non-reduced generators, reduced generators]
        options.R0 = zonotope([center(lastRtp),G(:,idxLength(gen:end)),...
            diag(sum(abs(G(:,idxLength(1:gen-1))),2))]);
        break
    end
    addErrorSprev = addErrorS;
end

options.eSi(options.i,1) = options.eS;

end

function options = errorAllocation(obj,options)
% separate two cases: no input, / input, time-varying inpu

if options.t == 0
    % init accumulating errors
    options.eS = 0;
    options.eSprev = 0;
    options.eP = 0;
    options.ePbox = 0;
    options.ePprev = 0;
end

if ~options.isInput
    % only homogeneous solution
    options.eHmax = options.error;
    options.ePmax = 0;
    options.eSmax = 0;
    
else %if ~isfield(options,'uTransVec')
    % inhomogeneous solution, but no uTransVec (no outer loop)
    options.eHmax = options.error * 0.30;
    options.ePmax = options.error * 0.30;
    options.eSmax = options.error * 0.40;
    
% else
%     % uTransVec -> outer loop
%     idxNextSwitch = find(options.uTransVecSwitch > options.t,1,'first');
%     if ~isempty(idxNextSwitch)
%         tFinalSwitch = options.uTransVecSwitch(idxNextSwitch);
%     else
%         tFinalSwitch = options.tFinal;
%     end
%     totalErrorSwitch = options.error / options.tFinal * ...
%         (tFinalSwitch - options.t);
%     % allocate homogeneous error
% %     options.eHmax = totalErrorSwitch * 0.30;
%     options.eHmax = options.error * 0.30;
%     
%     % allocate inhomogeneous error and reduction error
%     options.ePmax = options.error * 0.30;
%     options.eSmax = totalErrorSwitch * 0.40 ;
    
end

end

function [idxTimeStep, options] = checkTimeStep(obj,options)
% rA_postTSTT_checkTimeStep - check if current time step has been used
%    before and returns index for options.savedSets
% furthermore, if the time step is new and 0 \notin input set,
%    the necessary sets to compute the hom. and inhom. errors are computed

idxTimeStep = find(abs(options.usedTimeSteps - options.timeStep) < 1e-9,1); 

% time-varying input center: check if sets can be reused
% ...this is only possible if same input vectors are used
recompuTrans = false;
if ~isempty(idxTimeStep)
    if ~all(abs(options.savedSets{idxTimeStep,1}.uTrans - options.uTrans) < eps)
        recompuTrans = true;
    end
end

if isempty(idxTimeStep) || recompuTrans
    
    if ~recompuTrans
        % new time step
        idxTimeStep = length(options.savedSets) + 1;
        options.usedTimeSteps(idxTimeStep,1) = options.timeStep;
        options.savedSets{idxTimeStep,1}.timeStep = options.timeStep;
    end
    options.savedSets{idxTimeStep,1}.uTrans = options.uTrans;
    
    % initialize options.factor
    options.factor = 0;
    % used in inputTie and inputSolution
    for i=1:(options.taylorTerms+1)
        % compute initial state factor
        options.factor(i) = options.timeStep^(i)/factorial(i);    
    end
    options.savedSets{idxTimeStep}.factor = options.factor;
    
    if ~recompuTrans
        % compute error in exponential matrix
        % returns: obj.taylor.powers|error
        obj = exponential(obj, options);
        options.savedSets{idxTimeStep}.powers = obj.taylor.powers;
        options.savedSets{idxTimeStep}.E = obj.taylor.error;
    else
        obj.taylor.powers = options.savedSets{idxTimeStep}.powers;
        obj.taylor.error = options.savedSets{idxTimeStep}.E;
    end
    if options.isInput
        options.savedSets{idxTimeStep}.Etu = ...
            obj.taylor.error * options.timeStep * (obj.B * options.U);
    end

    if ~recompuTrans
        % compute time interval error (tie)
        % returns: obj.taylor.F
        obj = tie(obj,options);
        options.savedSets{idxTimeStep}.F = obj.taylor.F;
        options.savedSets{idxTimeStep}.FR0 = obj.taylor.F * options.R0;
    else
        % write to obj
        obj.taylor.F = options.savedSets{idxTimeStep}.F;
    end

    % obtain variables for inputTie
    options.savedSets{idxTimeStep}.inputCorr = zeros(obj.dim,1);
    if ~options.originContained
        % compute reachable set due to constant input vector:
        % returns: obj.taylor.inputF|inputCorr
        % ...calc. obj.taylor.V|RV|Rtrans|eAtInt|timeStep|eAt
        obj = inputTie(obj,options);
        % compute inputCorr 
        if isempty(obj.c); vTrans = obj.B*options.uTrans;
        else;              vTrans = obj.B*options.uTrans + obj.c; end
        obj.taylor.inputCorr = obj.taylor.inputF * zonotope(vTrans);
        options.savedSets{idxTimeStep}.inputCorr = obj.taylor.inputCorr;
    end
end

end

function maxTimeStep = maxTimeStepSize(obj,options)
% definition: max time step as time until input changes
% if no time-varying input vector, then return time until end of analysis

% non time-varying input center
maxTimeStep = options.tFinal - options.t;

% time-varying input center
if isfield(options,'uTransVec')
    idxNextSwitch = find(options.uTransVecSwitch > options.t,1,'first');
    if ~isempty(idxNextSwitch) % last step: idxNextSwitch = []
        maxTimeStep = options.uTransVecSwitch(idxNextSwitch) - options.t;
    end
end

end

function options = initTimeStep(obj,options)
% define shrinking factor and initialize time step size for loop
shrinkT = 1 / options.timeStepScaling; % (0, 1)
options.timeStep = maxTimeStepSize(obj,options);

U = obj.B * options.U;

% loop until bloatings due to Q * F * R0 and Q * E * t * u ...
%    below their respective allowed thresholds
while true
    
    % returns obj.taylor.error|powers|F: required for error computation
    try
        [obj,options] = expmtie_adaptive(obj,options);
    catch ME
        if strcmp(ME.identifier,'expmtie:notconverging')
            options.timeStep = options.timeStep * 0.5;
            continue;
        else
            rethrow(ME);
        end
    end
    
    % compute inhom. error
    if options.isInput
        % (re-)compute bound for QEtu
        fracRemTime = options.timeStep / (options.tFinal - options.t);
        ePadm = (options.ePmax - options.ePprev - options.eP) * fracRemTime;
        % compute error
        Etu = obj.taylor.error * options.timeStep * U;
        ePbox = options.ePbox + sum(abs(generators(Etu)),2);
        eP = norm(ePbox,2);
        % if inhom. error already to big, start loop anew
        if eP > ePadm
            options.timeStep = options.timeStep * shrinkT;
            continue;
        end
    end
    
    
	% compute hom. error
    FR0 = obj.taylor.F * options.R0;
    % obtain factors for inputSolution > inputTie (only if 0 not in input set)
    if ~options.originContained
        % compute reachable set due to input:
        % returns: obj.taylor.inputF|inputCorr
        % ...calc. obj.taylor.V|RV|Rtrans|eAtInt|timeStep|eAt in rA_initReach_new
        obj = inputTie(obj,options);
        if isempty(obj.c); vTrans = obj.B*options.uTrans;
        else;              vTrans = obj.B*options.uTrans + obj.c; end
        obj.taylor.inputCorr = obj.taylor.inputF * zonotope(vTrans);
        % compute induced bloating due to F * X0 + Ftilde * zon(uTrans)
        eH = norm(sum(abs(generators(FR0)),2),2) + ...
            norm(sum(abs(generators(obj.taylor.inputCorr)),2),2);
    else
        % compute induced bloating due to F * X0 (inputCorr = 0)
        eH = norm(sum(abs(generators(FR0)),2),2);
    end

    % check hom. error
    if eH > options.eHmax
        options.timeStep = options.timeStep * shrinkT;
        continue;
    end
    
    % thresholds satisfied -> exit
    break;
    
end

if options.i == 1
    idxTimeStep = 1;
elseif isfield(options,'uTransVec')    
    idxTimeStep = find(abs(options.usedTimeSteps - options.timeStep) < 1e-9,1); 
    if isempty(idxTimeStep)
        idxTimeStep = length(options.usedTimeSteps)+1;
    end
end

% saved sets for look-up instead of recalculation (time save)
options.usedTimeSteps(idxTimeStep) = options.timeStep;
options.savedSets{idxTimeStep}.timeStep = options.timeStep;
options.savedSets{idxTimeStep}.uTrans = options.uTrans;
options.savedSets{idxTimeStep}.powers = obj.taylor.powers;
options.savedSets{idxTimeStep}.E = obj.taylor.error;
if options.isInput
    options.savedSets{idxTimeStep}.Etu = Etu;
end
options.savedSets{idxTimeStep}.F = obj.taylor.F;
options.savedSets{idxTimeStep}.FR0 = FR0;
options.savedSets{idxTimeStep}.factor = options.factor;
if ~options.originContained
    options.savedSets{idxTimeStep}.inputCorr = obj.taylor.inputCorr;
else
    options.savedSets{idxTimeStep}.inputCorr = zeros(obj.dim,1);
end

% exponential matrix for given time step size
eAt = expm(obj.A*options.timeStep);

% save current propagation matrix in same struct
obj.taylor.eAt = eAt;

% initialize propagation matrices for analysis
options.P = eye(obj.dim);
options.Q = obj.taylor.eAt;

options.eHi(options.i,1) = eH;
% accumulate bloating
if options.isInput
    options.ePbox = ePbox;
    options.eP = eP;
    options.ePi(options.i,1) = options.ePprev + eP;
end

end

function options = rA_initReach_new(obj,options)
% Inputs:
%    obj      - linearSys object
%    options  - options struct, containing
%                    .R0: initial set
%                    .timeStep: time step size for current iteration
%                    .taylorTerms: taylor terms (fixed)
%
% Outputs:
%    options - options struct, containing
%                    .Rhom: homogeneous time-interval solution
%                    .Rhom_tp: homogeneous time-point solution
%                    .Raux: inhomogeneous solution due to U
%                    .Rtrans: inhomogeneous solution due to uTrans
%                    .Rinhom: full inhomogeneous solution (only first step)

% compute reachable set due to input:
% returns: obj.taylor.V|RV|Rtrans|inputF|inputCorr|eAtInt
obj = inputSolution(obj,options);

% inhomogeneous solution due to U
RV = zonotope(zeros(obj.dim,1));
if any(any(options.U.Z))
    gensRV = generators(obj.taylor.RV);
    RV = zonotope([zeros(obj.dim,1), gensRV(:,any(gensRV,1))]);
end
% inhomogeneous solution due to u
Rtrans = zonotope(zeros(obj.dim,1));
if any(options.uTrans)
    gensRtrans = generators(obj.taylor.Rtrans);
    Rtrans = zonotope([center(obj.taylor.Rtrans), ...
        gensRtrans(:,any(gensRtrans,1))]);
end

% first time step homogeneous solution
Rhom_tp = obj.taylor.eAt*options.R0 + Rtrans;
Rhom = enclose(options.R0,Rhom_tp) + obj.taylor.F*options.R0 + ...
    obj.taylor.inputCorr;


% save homogeneous solution
options.Rhom    = Rhom;
options.Rhom_tp = Rhom_tp;
% save inhomogeneous solution
options.Raux    = RV;
options.Rtrans  = Rtrans;

% first step (per uTrans): init handling of inhom. solution
if options.i == 1 || (isfield(options,'uTransVec') && ...
        any(abs(options.t - options.uTransVecSwitch) < 1e-9) )
    options.lastRtrans = options.Rtrans;
    options.Rinhom = RV;

    % prepare data for reduction of inhom. solution in main loop
    if options.isInput
        % generators of Rinhom
        options.RinhomG = generators(options.Rinhom);
        options.RinhomGlength = vecnorm(options.RinhomG,2);
        options.RinhomRed = 0;
        % container for generators converted to interval (always diag matrix)
        options.RinhomConvInt = zeros(obj.dim);
    end
end

end

function options = postTimeStep(obj,options)

maxTimeStep = maxTimeStepSize(obj,options);
if options.timeStep * options.timeStepScaling > maxTimeStep
    options.timeStep = maxTimeStep;
else
    % try bigger time step
    options.timeStep = options.timeStep * options.timeStepScaling;
end


% loop until bloatings due to Q * F * R0 and Q * E * t * u ...
%    below their respective allowed thresholds
while true
    
    % recompute bound for QEtu
    if options.isInput
        fracRemTime = options.timeStep / (options.tFinal - options.t);
        % allowed error for this step: allowed - curr divided by forseen time
        ePadm = (options.ePmax - options.ePprev - options.eP) * fracRemTime;
        % ... if no particular solution, then ePadm same as in init
    end
    % check if time step already used
    [idxTimeStep, options] = checkTimeStep(obj,options);
    % ... if not, required sets computed
    
    % compute hom. error
    eH = norm(sum(abs(generators(options.Q * ...
        options.savedSets{idxTimeStep}.FR0)),2),2);
    if ~options.originContained
        eH = eH + norm(sum(abs(generators(options.Q * ...
            options.savedSets{idxTimeStep}.inputCorr)),2),2);
    end
    
    % check hom. error
    if eH > options.eHmax
        options.timeStep = options.timeStep / options.timeStepScaling;
        continue;
    end
    
    % compute inhom. error
    if options.isInput
        Etu = options.Q * options.savedSets{idxTimeStep}.Etu;
        ePbox = options.ePbox + sum(abs(generators(Etu)),2);
        eP = norm(ePbox,2);
        % check inhom. error
        if eP > ePadm
            options.timeStep = options.timeStep / options.timeStepScaling;
            continue;
        end
    end
    
    % both errors satisfied
    break;

end


% for inputSolution
if ~options.originContained
    options.factor = options.savedSets{idxTimeStep}.factor;
    obj.taylor.inputCorr = options.savedSets{idxTimeStep}.inputCorr;
end

obj.taylor.powers = options.savedSets{idxTimeStep}.powers;
obj.taylor.error = options.savedSets{idxTimeStep}.E;
obj.taylor.F = options.savedSets{idxTimeStep}.F;
% obj.taylor is correct with the exception of:
%    obj.taylor.V|RV|Rtrans|eAtInt|timeStep|eAt
% ...currently, obj.taylor.inputF|inputCorr are also calculated later
obj.taylor.eAt = expm(obj.A*options.timeStep);


options.eHi(options.i,1) = eH;
% accumulate bloating
if options.isInput
    options.ePbox = ePbox;
    options.eP = eP;
    options.ePi(options.i,1) = options.ePprev + eP;
end


end

function options = propInhom(obj,options)
% rA_propInhom - propagate inhomogeneous solution,
%    including reduction of zonotope order until limit given by
%    allowed bloating reached
% formula given by:
% P^R([0, t_i + Delta t_i]) = P^R([0, t_i])
%    + e^At_i * P^R_U(Delta t_i) + e^At_i-1 * P^R_u(Delta t_i-1)

% calculate linearly growing boundary for allowed bloating error
fracRemTime = options.timeStep / (options.tFinal - options.t);
eSadm = (options.eSmax - options.eSprev - options.eS) * fracRemTime;

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
        if options.eS + addErrorS > eSadm
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
        if options.eS + addErrorS > eSadm || gen == numGensInhom
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

options.eSi(options.i,1) = options.eSprev + options.eS;

end

%------------- END OF CODE --------------