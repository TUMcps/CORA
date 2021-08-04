function pZ = restoreId(pZ,id)
% restoreId - adds all ids in "id" to pZ (if not already present)
%
% Syntax:  
%    res = restoreId(pZ,id)
%
% Inputs:
%    pZ  - polyZonotope object
%    id  - ids to be restored
%
% Outputs:
%    pZ - resulting polynomial Zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: replaceId

% Author:       Victor Gassmann
% Written:      13-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
ind = ismember(id,pZ.id);
n_toRestore = sum(~ind);
ne = size(pZ.expMat,2);
pZ.expMat = [pZ.expMat;zeros(n_toRestore,ne)];
pZ.id = [pZ.id;id(~ind)];
%------------- END OF CODE --------------