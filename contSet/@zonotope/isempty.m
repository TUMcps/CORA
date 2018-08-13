function res = isempty(obj)
% is_empty - returns 1 if a zonotope is empty and 0 otherwise
%
% Syntax:  
%    res = is_empty(obj)
%
% Inputs:
%    obj - zonotope object
%
% Outputs:
%   res - result in {0,1}
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      21-August-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = isempty(obj.Z);

%------------- END OF CODE --------------