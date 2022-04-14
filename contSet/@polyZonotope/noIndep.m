function PZ = noIndep(pZ)
% noIndep - removes Grest from pZ
%
% Syntax:  
%    PZ = noIndep(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    PZ - polyZonotope object
%
% Example: 
%    pZ1 = polyZonotope([2;1],[1 0; -2 1],[1; 0],[0 2; 1 0]);
%    noIndep(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Victor Gassmann
% Written:       21-January-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
PZ = polyZonotope(pZ.c,pZ.G,zeros(length(pZ.c),0),pZ.expMat,pZ.id);
%------------- END OF CODE --------------