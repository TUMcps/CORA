function val = distance(obj,S)
% distance - computes the distance from a conHyperplane object to the set S
%
% Syntax:  
%    res = distance(obj,S)
%
% Inputs:
%    obj - conHyperplane object
%    S   - set
%
% Outputs:
%    val - distance between obj and S
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      08-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if isa(S,'double')
    val = distancePoint(obj,S);
else
    error('Second argument type not supported');
end
%------------- END OF CODE --------------