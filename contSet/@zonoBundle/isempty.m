function res = isempty(obj)
% isempty - checks if zonoBundle is empty or not
%
% Syntax:  
%    res = isempty(obj)
%
% Inputs:
%    obj - zonoBundle object
%
% Outputs:
%    res - boolean whether obj is empty or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger, Niklas Kochdumper
% Written:      16-Sep-2019
% Last update:  21-Nov-2019 (NK, call conZonotope function)
% Last revision:---

%------------- BEGIN CODE --------------

    % check if any of the single zonotopes is empty
    for i = 1:length(obj.Z)
       if isempty(obj.Z{i})
           res = 1;
           return;
       end
    end
    
    % check if the intersection of the zonotopes is empty
    res = isempty(conZonotope(obj));

%------------- END OF CODE --------------