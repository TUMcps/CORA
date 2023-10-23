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
%    Z = zonotope([0;0],[1 3 2 -2; 3 -2 0 1]);
%    Glength = generatorLength(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       19-July-2010
% Last update:   14-March-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Glength = vecnorm(Z.G);

% ------------------------------ END OF CODE ------------------------------
