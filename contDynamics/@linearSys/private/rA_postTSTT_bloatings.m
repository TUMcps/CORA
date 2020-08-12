function [inducedBloatingF, addBloatingInhom] = rA_postTSTT_bloatings(obj, options, idxTimeStep)
% rA_postTSTT_bloatings - calculate bloatings due to
%   Q * F * X0 + Q * inputCorr
%   Q * E * timeStep * u ... with u = uTrans + U
% using interval overapproximation
%
% Syntax:  
%    [inducedBloatingF, addBloatingInhom] = ...
%           rA_postTSTT_bloatings(obj, options, idxTimeStep)
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

if ~options.originContained
    inducedBloatingF = ...
        norm(sum(abs(generators(options.Q * ...
            options.savedSets{idxTimeStep}.FR0{options.taylorTerms})),2),2) + ...
        norm(sum(abs(generators(options.Q * ...
            options.savedSets{idxTimeStep}.inputCorr{options.taylorTerms})),2),2);
else
    inducedBloatingF = ...
        norm(sum(abs(generators(options.Q * ...
            options.savedSets{idxTimeStep}.FR0{options.taylorTerms})),2),2);
end

if options.isInput
    addBloatingInhom = ...
        norm(sum(abs(generators(options.Q * ...
        options.savedSets{idxTimeStep}.Etu{options.taylorTerms})),2),2);
else
    addBloatingInhom = norm(sum(abs(generators(options.Q * ...
        (options.savedSets{idxTimeStep}.E{options.taylorTerms} * options.R0))),2),2);
end


end

%------------- END OF CODE --------------
