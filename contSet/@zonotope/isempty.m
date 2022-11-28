function res = isempty(Z)
% isempty - checks if a zonotope is the empty set
%
% Syntax:  
%    res = isempty(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z1 = zonotope();
%    Z2 = zonotope([0;0]);
%    Z3 = zonotope([1;0],[1 -2 1; 0 1 -2]);
%
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

res = isempty(Z.Z) && isempty(Z.halfspace);

%------------- END OF CODE --------------