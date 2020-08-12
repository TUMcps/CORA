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
    obj = rA_inputTie(obj,options);
    options.savedSets{idxTimeStep}.inputCorr{options.taylorTerms,1} = ...
        obj.taylor.inputCorr;
else
    options.savedSets{idxTimeStep}.inputCorr{options.taylorTerms,1} = ...
        zeros(obj.dim,1);
end
% -------------------------------------------------------------------------



end

%------------- END OF CODE --------------
