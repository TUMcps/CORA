function [pZ_x,pZ_nx] = onlyId(pZ,idx)
% onlyId - splits pZ into part only dependent on 'id' and rest
%
% Syntax:  
%    res = onlyId(pZ,id)
%
% Inputs:
%    pZ - polyZonotope object
%    id- id
%
% Outputs:
%    res - 1 if set is all zero, 0 otherwise
%
% Example: 
%    pZ = polyZonotope([0;0],[1,1],zeros(2,0),[1,2;3,4;5,6]);
%    res = isZero(pZ);    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope, isPolytope

% Author:       Victor Gassmann
% Written:      12-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
ind_x = ismember(pZ.id,idx);
n = length(pZ.c);
if sum(ind_x)<length(idx)
    error('id contains ids not contained in pZ.id');
end
ind_mask = sum(pZ.expMat(~ind_x,:),1)==0;
G_x = pZ.G(:,ind_mask);
eM_x = pZ.expMat(ind_x,ind_mask);
G_nx = pZ.G(:,~ind_mask);
eM_nx = pZ.expMat(:,~ind_mask);
pZ_x = polyZonotope(pZ.c,G_x,zeros(n,0),eM_x,pZ.id(ind_x));
pZ_nx = polyZonotope(zeros(n,1),G_nx,pZ.Grest,eM_nx,pZ.id);
%------------- END CODE --------------