function c = center(Z)
% center - returns the center of a zonotope
%
% Syntax:  
%    c = center(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    c - center of the zonotope Z
%
% Example:
%    Z = zonotope([1;0],[1 0; 0 1]);
%    c = center(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      30-September-2006 
% Last update:  22-March-2007
%               14-March-2021 (MW, empty set)
% Last revision:---

%------------- BEGIN CODE --------------

try
    c = Z.Z(:,1);
catch ME
    if isempty(Z)
        c = [];
    else
        rethrow(ME);
    end
end
    

%------------- END OF CODE --------------