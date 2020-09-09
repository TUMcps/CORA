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
%    obj.taylor.V|RV|Rtrans|eAtInt|timeStep|eAt (all in syn_post)
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

%------------- END OF CODE --------------
