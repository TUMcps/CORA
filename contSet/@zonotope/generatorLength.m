function [genLen] = generatorLength(Z)
% generatorLength - Returns the lengths of the generators
%
% Syntax:  
%    [genLen] = generatorLength(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    genLen - vector of generator length
%
% Example: ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      19-July-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

genLen=vnorm(Z.Z(:,2:end),1);

%------------- END OF CODE --------------