function pZ = deleteZeros(pZ)
% deleteZeros - deletes all generators of length 0
%
% Syntax:
%    pZ = deleteZeros(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example:
%    pZ = polyZonotope([1;3],[1 2 0;1 -2 0],[1 0 2;0 0 -1],[1 2;1 1;2 1]);
%    pZ = deleteZeros(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/deleteZeros

% Author:       Mark Wetzlinger
% Written:      20-April-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% indices with non-zero generators
idxD = any(pZ.G,1);
idxI = any(pZ.Grest,1);

% if all non-zero, skip
if ~( all(idxD) && all(idxI) ) 
    % delete zero generators
    pZ.G = pZ.G(:,idxD);
    pZ.expMat = pZ.expMat(:,idxD);
    % something to do with pZ.id?
    pZ.Grest = pZ.Grest(:,idxI);
else
    return
end



%------------- END OF CODE --------------