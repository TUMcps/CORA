function res = isempty(obj)
% isempty - returns true if a zonotope is empty and false otherwise
%
% Syntax:  
%    res = isempty(obj)
%
% Inputs:
%    obj - zonotope object
%
% Outputs:
%    res - result in {0,1}
%
% Example: 
%    Z1 = zonotope([]);
%    Z2 = zonotope([0;0]);
%    Z3 = zonotope([1;0],[1 -2 1; 0 1 -2]);
%    isempty(Z1); % true
%    isempty(Z2); % false
%    isempty(Z3); % false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      21-August-2015
% Last update:  01-May-2020 (MW, add other properties of obj)
% Last revision:---

%------------- BEGIN CODE --------------

res = isempty(obj.Z) && isempty(obj.halfspace);

%------------- END OF CODE --------------