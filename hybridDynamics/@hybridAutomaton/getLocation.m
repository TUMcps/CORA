function [location] = getLocation(obj, locIndex)
% getLocation - Get the location object from flat automaton obj at index locIndex.
%
% Syntax:  
%    location = getLocation(obj, locIndex)
%
% Input:
%     locIndex - Index for a location (numeral or 1x1 cell)
% Author:       Johann Schöpfer
% Written:      24-May-2016  
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------
if isnumeric(locIndex)
    location = obj.location{locIndex};
elseif iscell(locIndex) && length(locIndex) == 1 && isnumeric(locIndex{1})
    location = obj.location{locIndex{1}};
else
    error("Location index for a flat automaton must be numeric or 1x1 cell");
end

