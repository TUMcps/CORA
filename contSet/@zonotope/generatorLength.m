function Glength = generatorLength(Z)
% generatorLength - returns the lengths of the generators
%
% Syntax:  
%    Glength = generatorLength(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Glength - vector of generator length
%
% Example:
%    Z = zonotope([0;0],rand(2,10));
%    Glength = generatorLength(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      19-July-2010
% Last update:  14-March-2019
% Last revision:---

%------------- BEGIN CODE --------------

Glength = vecnorm(Z.Z(:,2:end));

%------------- END OF CODE --------------