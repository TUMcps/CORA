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
        obj = rA_inputTie(obj,options);
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
        obj = rA_inputTie(obj,options);
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


%------------- END OF CODE --------------