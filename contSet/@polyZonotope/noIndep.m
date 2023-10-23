function pZ = noIndep(pZ)
% noIndep - removes independent generator matrix from polynomial zonotope
%
% Syntax:
%    pZ = noIndep(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    pZ = polyZonotope([2;1],[1 0; -2 1],[1; 0],[0 2; 1 0]);
%    noIndep(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       21-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

pZ.GI = zeros(dim(pZ),0);

% ------------------------------ END OF CODE ------------------------------
