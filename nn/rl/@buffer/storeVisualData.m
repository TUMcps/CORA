function obj = storeVisualData(obj,reachSet,episodeNum)
% storeVisualData - stor visual data for rendering of environment
%
% Syntax:
%   obj = storeVisualData(obj,reachSet,episodeNum)
%
% Inputs:
%   reachSet - reachable set for environment step
%   episodeNum - number of the current episode
%
% Outputs:
%   obj - updated buffer
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: buffer

% Authors:       Manuel Wendl
% Written:       03-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({ ...
    {reachSet, 'att', {'struct','reachSet'}}, ...
    {episodeNum, 'att', 'numeric'}
    })

if isempty(obj.visualisationData.episodeNum)
    obj.visualisationData.episodeNum = episodeNum;
    obj.visualisationData.reachSet = {reachSet};
elseif episodeNum > obj.visualisationData.episodeNum(end)
    obj.visualisationData.episodeNum = [obj.visualisationData.episodeNum, episodeNum];
    obj.visualisationData.reachSet{end+1} = reachSet;
else
    if isa(reachSet,'struct')
        obj.visualisationData.reachSet{end}.t = [obj.visualisationData.reachSet{end}.t;reachSet.t];
        obj.visualisationData.reachSet{end}.x = [obj.visualisationData.reachSet{end}.x;reachSet.x];
    else
        obj.visualisationData.reachSet{end} = add(obj.visualisationData.reachSet{end},reachSet);
    end
end
end

% ------------------------------ END OF CODE ------------------------------
