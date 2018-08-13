function vol = volume(matI)
% volume - computes the volume of an interval matrix by computing the 
% volume of the corresponding interval hull
%
% Syntax:  
%    vol = volume(matI)
%
% Inputs:
%    matI - interval matrix
%
% Outputs:
%    vol - volume
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      24-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%conversion to a zonotope
IH = intervalhull(matI);

%compute volume of the zonotope
vol = volume(IH);

%------------- END OF CODE --------------