function res = isempty(zB)
% isempty - checks if a zonotope bundle is the empty set
%
% Syntax:  
%    res = isempty(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    res - true/false
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
for i = 1:length(zB.Z)
   if isempty(zB.Z{i})
       res = true;
       return;
   end
end

% check if the intersection of the zonotopes is empty
res = isempty(conZonotope(zB));

%------------- END OF CODE --------------