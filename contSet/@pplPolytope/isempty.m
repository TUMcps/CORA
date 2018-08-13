function res = isempty(obj)
% is_empty - returns 1 if a pplPolytope is empty and 0 otherwise
%
% Syntax:  
%    res = is_empty(obj)
%
% Inputs:
%    obj - pplPolytope object
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
% Written:      20-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute vertices
res = is_empty(-obj.C,obj.d);


%------------- END OF CODE --------------