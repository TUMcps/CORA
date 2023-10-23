function pZ = rmId(pZ,id)
% rmId - removes identifiers from polynomial zonotope
%
% Syntax:
%    pZ = rmId(pZ,id)
%
% Inputs:
%    pZ - polyZonotope object
%    id - identifiers that should to be removed
%
% Outputs:
%    pZ - polyZonotope object
%
% Example:
%    pZ = polyZonotope([1;2],[1,2;3,4],[],[1,2;2,1],[1;2]);
%    % visualize
%    print(pZ,'Ids',{1,2},'Vars',{'x','y'});
%    print(rmId(pZ,1),'Ids',{2},'Vars',{'y'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Victor Gassmann
% Written:       16-February-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[~,ii_rm] = intersect(pZ.id,id,'stable');
pZ.E(ii_rm,:) = [];
pZ.id(ii_rm) = [];

% remove any possibly-generated same-exponent vector entries
[pZ.E,pZ.G] = removeRedundantExponents(pZ.E,pZ.G);

% ------------------------------ END OF CODE ------------------------------
