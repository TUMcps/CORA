function [activeGuards] = hyperplaneIntersectCheck(obj,R)
% hyperplaneIntersectCheck - checks if zonotope intersects a hyperplane and
% returns the guard number
%
% Syntax:  
%    [activeGuards] = hyperplaneIntersectCheck(obj,R)
%
% Inputs:
%    obj - location object
%    R - reachable set
%
% Outputs:
%    activeGuards - guards that activate transitions
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      17-October-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%initialize values
activeGuards=[]; 

%check if guard sets exist
if ~isempty(obj.transition)
    for iTransition=1:length(obj.transition)
        guard = get(obj.transition{iTransition},'guard');
        if zonoIntersect(guard, R);
            activeGuards(end+1) = iTransition; 
        end
    end
end

%------------- END OF CODE --------------