function res = isempty(obj)
% is_empty - returns 1 if there is no vertex and 0 otherwise
%
% Syntax:  
%    res = is_empty(obj)
%
% Inputs:
%    obj - vertices object
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
% Written:      20-April-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = isempty(obj.V);

%------------- END OF CODE --------------