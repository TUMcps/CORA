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
    % used in rA_inputTie and inputSolution
    for i=1:(options.maxTaylorTerms+1)
        % compute initial state factor
        options.factor(i) = options.timeStep^(i)/factorial(i);    
    end
    options.savedSets{idxTimeStep}.factor = options.factor;
end



end

%------------- END OF CODE --------------

