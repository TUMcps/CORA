function r = enclosingRadius(obj)
% enclosingRadius - Computes radius of enclosing hyperball of an interval 
%
% Syntax:  
%    r = enclosingRadius(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    r - radius of enclosing hyperball
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
% Written:      22-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compte interval radius
intervalRadius = rad(obj);

%compute radius
r = sqrt(sum(intervalRadius.^2));

%------------- END OF CODE --------------